#ifdef WITH_NLW2

#include <unordered_set>

#include "mp/nl-solver.h"
#include "mp/nl-feeder.h"
#include "mp/sol-handler.h"
#include "mp/nl-opcodes.h"
#include "mp/nl-solver.hpp"

extern "C" {
#include "mpl.h"

struct csa *mathprog_nlw2_build_prob(mpl_tran *tran, int relax);
int mathprog_nlw2_postsolve(mpl_tran *tran, struct csa *csa, int relax);
int mathprog_nlw2_run_solver(mpl_tran *tran, struct csa *csa, const char *solver, const char *options);
void mathprog_nlw2_cleanup(struct csa *csa);
}

struct csa
{     /* common storage area */
      int m, n;
      char *prob_name;
      char **var_names, **row_names;
      double *col_size;
      size_t nnz;
      int solve_result;
      double *x, *y, *var_status, *row_status;
      int m_range, m_eqn, n_bin, n_int, nnz_obj;
      int *perm, *inv_perm;
      int obj_no;
      int *obj_row, *con_row;
      int *ind;
      double *val;
      int m_obj, m_con;
      int n_nl_obj, n_nl_obj_i;
      int n_nl_con, n_nl_con_i;
      int n_nl_both, n_nl_both_i;
      int m_nl, m_nl_obj;
      int n_def;
      int *defvar_order;
};

static void find_def_var_used(NL_EXPR *nl, std::unordered_set<int> &used)
{
   if (nl == NULL) return;
   if (nl->op == O_LINFORM) return;
   if (is_nary_nl_expr_op(nl->op))
   {  for (int i = 0; i < nl->u.vec.n; i++)
         find_def_var_used(nl->u.vec.data[i], used);
      return;
   }
   if (nl->op != O_DEFVAR)
   {  find_def_var_used(nl->u.arg.x, used);
      find_def_var_used(nl->u.arg.y, used);
      find_def_var_used(nl->u.arg.z, used);
      return;
   }
   used.insert(nl->u.defvar->j-1);
   return;
}


struct csa *mathprog_nlw2_build_prob(mpl_tran *tran, int relax)
{     /* build LP/MIP problem instance from the model */
      struct csa *csa = new struct csa;
      int m, n, i, j, t, ct, kind, type, len, *ind, nnz, n_def;
      double lb, ub, *val;
      char *name;
      int *cnt_idx, *int_idx, *bin_idx, *nl_both_idx, *nl_both_int_idx;
      int *nl_obj_idx, *nl_obj_int_idx, *nl_con_idx, *nl_con_int_idx;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mpl_build_prob: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      n_def = mpl_get_num_defs(tran);
      csa->ind = ind = (int *)xcalloc(1+n, sizeof(int));
      cnt_idx = (int *)xcalloc(n, sizeof(int));
      int_idx = (int *)xcalloc(n, sizeof(int));
      bin_idx = (int *)xcalloc(n, sizeof(int));
      nl_obj_idx = (int *)xcalloc(n, sizeof(int));
      nl_obj_int_idx = (int *)xcalloc(n, sizeof(int));
      nl_con_idx = (int *)xcalloc(n, sizeof(int));
      nl_con_int_idx = (int *)xcalloc(n, sizeof(int));
      nl_both_idx = (int *)xcalloc(n, sizeof(int));
      nl_both_int_idx = (int *)xcalloc(n, sizeof(int));
      csa->val = val = (double *)xcalloc(1+n, sizeof(double));
      csa->m = m;
      csa->n = n;
      csa->n_def = n_def;
      csa->var_names = (char **)xcalloc(n + n_def, sizeof(char *));
      csa->row_names = (char **)xcalloc(m, sizeof(char *));
      csa->col_size = (double *)xcalloc(n, sizeof(double));
      csa->solve_result = -1;
      csa->x = NULL;
      csa->y = NULL;
      csa->var_status = NULL;
      csa->row_status = NULL;
      csa->m_range = 0;
      csa->m_eqn = 0;
      csa->n_bin = 0;
      csa->n_int = 0;
      csa->nnz_obj = 0;
      csa->perm = (int *)xcalloc(n, sizeof(int));
      csa->inv_perm = (int *)xcalloc(n, sizeof(int));
      csa->obj_row = (int *)xcalloc(m, sizeof(int));
      csa->con_row = (int *)xcalloc(m, sizeof(int));
      csa->obj_no = -1;
      csa->m_obj = 0;
      csa->m_con = 0;
      csa->n_nl_obj = 0;
      csa->n_nl_obj_i = 0;
      csa->n_nl_con = 0;
      csa->n_nl_con_i = 0;
      csa->n_nl_both = 0;
      csa->n_nl_both_i = 0;
      csa->m_nl = 0;
      csa->m_nl_obj = 0;
      csa->defvar_order = NULL;

      ct = 0;
      for (j = 0; j < n; j++)
      {  name = mpl_get_col_name(tran, j+1);
         csa->var_names[j] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->var_names[j], name);
         kind = mpl_get_col_kind(tran, j+1);
         csa->col_size[j] = 0;
         /* split variables per kind for permutation array */
         switch(mpl_get_col_nl_flags(tran, j+1))
         {  case MPL_LIN:
               if (kind == MPL_NUM || relax)
                  cnt_idx[ct++] = j;
               else if (kind == MPL_BIN)
                  bin_idx[csa->n_bin++] = j;
               else
                  int_idx[csa->n_int++] = j;
               break;
            case MPL_NLO:
               if (kind == MPL_NUM || relax)
                  nl_obj_idx[csa->n_nl_obj++] = j;
               else
                  nl_obj_int_idx[csa->n_nl_obj_i++] = j;
               break;
            case MPL_NLC:
               if (kind == MPL_NUM || relax)
                  nl_con_idx[csa->n_nl_con++] = j;
               else
                  nl_con_int_idx[csa->n_nl_con_i++] = j;
               break;
            case MPL_NLO | MPL_NLC:
               if (kind == MPL_NUM || relax)
                  nl_both_idx[csa->n_nl_both++] = j;
               else
                  nl_both_int_idx[csa->n_nl_both_i++] = j;
               break;
            default:
               assert(kind != kind);
         }
      }
      /* the varialbles should be ordered nl_both, nl_both_int, nl_cont,
         nl_cont_int, nl_obj, nl_obj_int, cnt, bin, int */
      memcpy(csa->perm, nl_both_idx, csa->n_nl_both * sizeof(int));
      t = csa->n_nl_both;
      memcpy(csa->perm + t, nl_both_int_idx, csa->n_nl_both_i * sizeof(int));
      t += csa->n_nl_both_i;
      memcpy(csa->perm + t, nl_con_idx, csa->n_nl_con * sizeof(int));
      t += csa->n_nl_con;
      memcpy(csa->perm + t, nl_con_int_idx, csa->n_nl_con_i * sizeof(int));
      t += csa->n_nl_con_i;
      memcpy(csa->perm + t, nl_obj_idx, csa->n_nl_obj * sizeof(int));
      t += csa->n_nl_obj;
      memcpy(csa->perm + t, nl_obj_int_idx, csa->n_nl_obj_i * sizeof(int));
      t += csa->n_nl_obj_i;
      memcpy(csa->perm + t, cnt_idx, ct * sizeof(int));
      t += ct;
      memcpy(csa->perm + t, bin_idx, csa->n_bin * sizeof(int));
      t += csa->n_bin;
      memcpy(csa->perm + t, int_idx, csa->n_int * sizeof(int));
      t += csa->n_int;
      assert(t == csa->n);
      /* initialize inverse permutation */
      for (j = 0; j < n; j++)
         csa->inv_perm[csa->perm[j]] = j;

      for (j = 0; j < n_def; j++)
      {  name = mpl_get_col_name(tran, n+j+1);
         csa->var_names[n+j] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->var_names[n+j], name);
      }

      nnz = 0;
      cnt_idx = (int *)xrealloc(cnt_idx, m, sizeof(int));
      for (i = 0; i < m; i++)
      {  name = mpl_get_row_name(tran, i+1);
         csa->row_names[i] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->row_names[i], name);
         kind = mpl_get_row_kind(tran, i+1);
         if (kind == MPL_MIN || kind == MPL_MAX)
         {  len = mpl_get_mat_row(tran, i+1, NULL, NULL);
            csa->nnz_obj += len;
            csa->obj_row[csa->m_obj] = i;
            csa->m_obj++;
            if (mpl_get_row_nl(tran, i+1))
               csa->m_nl_obj++;
         }
         else
         {  type = mpl_get_row_bnds(tran, i+1, NULL, NULL);
            if (type == MPL_FX)
               csa->m_eqn++;
            else if (type == MPL_DB)
               csa->m_range++;
            if (mpl_get_row_nl(tran, i+1))
               csa->con_row[csa->m_nl++] = i;
            else
               cnt_idx[csa->m_con++] = i;
            len = mpl_get_mat_row(tran, i+1, ind, NULL);
            nnz += len;
            for (t = 1; t <= len; t++)
               if (ind[t] <= n)
                  csa->col_size[ind[t]-1]++;
         }
      }
      memcpy(csa->con_row + csa->m_nl, cnt_idx, csa->m_con * sizeof(int));
      csa->m_con += csa->m_nl;
      assert(csa->m_con + csa->m_obj == csa->m);
      csa->nnz = nnz;

      name = mpl_get_prob_name(tran);
      csa->prob_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
      strcpy(csa->prob_name, name);

      if (n_def > 0)
      {  /* order defined variables to avoid forward references */
         csa->defvar_order = (int *)xcalloc(n_def, sizeof(int));
         std::vector<std::pair<int, std::unordered_set<int>>> defvar_usage;
         for (j = 0; j < n_def; j++)
         {  defvar_usage.emplace_back(n+j, std::unordered_set<int>{});
            find_def_var_used(tran->col[n+j+1]->def->nl,
               defvar_usage.back().second);
         }
         ct=0;
         while (!defvar_usage.empty())
         {  bool flag = false;
            for (auto it = defvar_usage.begin(); it != defvar_usage.end();)
            {  if (it->second.empty())
               {  int val = it->first;
                  it = defvar_usage.erase(it);
                  csa->defvar_order[ct++] = val;
                  for (auto &t : defvar_usage)
                     t.second.erase(val);
                  flag = true;
               }
               else
                  it++;
            }
            if (!flag) xerror("can't resolve defined variables order\n");
         }
         assert(ct == n_def);
      }

      /* free working arrays */
      xfree(cnt_idx);
      xfree(bin_idx);
      xfree(int_idx);
      xfree(nl_obj_idx);
      xfree(nl_obj_int_idx);
      xfree(nl_con_idx);
      xfree(nl_con_int_idx);
      xfree(nl_both_idx);
      xfree(nl_both_int_idx);

      return csa;
}

class Feeder : public mp::NLFeeder<struct csa*, void*>
{
public:
   Feeder(mpl_tran *t, struct csa *c) : tran(t), csa(c)
   {  hdr.num_vars = csa->n;
      hdr.num_algebraic_cons = csa->m_con;
      hdr.num_objs = csa->m_obj;
      hdr.num_ranges = csa->m_range;
      hdr.num_eqns = csa->m_eqn;
      hdr.num_linear_binary_vars = csa->n_bin;   //binary variables
      hdr.num_linear_integer_vars = csa->n_int;  //integer variables
      hdr.prob_name = csa->prob_name;
      hdr.format = mp::NLHeader::TEXT;
      hdr.flags = WANT_OUTPUT_SUFFIXES;
      hdr.num_con_nonzeros = csa->nnz;
      hdr.num_obj_nonzeros = csa->nnz_obj;

      // apparently those don't need to be set
      hdr.max_con_name_len = 0;
      hdr.max_var_name_len = 0;

      // support for nonlinear expressions
      hdr.num_nl_cons = csa->m_nl;
      hdr.num_nl_objs = csa->m_nl_obj;
      hdr.num_nl_vars_in_cons = csa->n_nl_both + csa->n_nl_con
         + csa->n_nl_both_i + csa->n_nl_con_i;
      // despite what the docs say, num_nl_vars_in_objs includes all
      // constraint vars
      hdr.num_nl_vars_in_objs = hdr.num_nl_vars_in_cons + csa->n_nl_obj
         + csa->n_nl_obj_i;
      hdr.num_nl_vars_in_both = csa->n_nl_both + csa->n_nl_both_i;
      hdr.num_nl_integer_vars_in_both = csa->n_nl_both_i;
      hdr.num_nl_integer_vars_in_cons = csa->n_nl_con_i;
      hdr.num_nl_integer_vars_in_objs = csa->n_nl_obj_i;

      // the following are not (yet?) supported by MathProg
      hdr.num_logical_cons = 0;
      hdr.num_compl_conds = 0;
      hdr.num_nl_compl_conds = 0;
      hdr.num_compl_dbl_ineqs = 0;
      hdr.num_compl_vars_with_nz_lb = 0;
      hdr.num_nl_net_cons = 0;
      hdr.num_linear_net_cons = 0;
      hdr.num_linear_net_vars = 0;
      hdr.num_funcs = 0;
      hdr.num_common_exprs_in_both = csa->n_def;
      hdr.num_common_exprs_in_cons = 0;
      hdr.num_common_exprs_in_objs = 0;
      hdr.num_common_exprs_in_single_cons = 0;
      hdr.num_common_exprs_in_single_objs = 0;
   }

   mp::NLHeader Header() { return hdr; }

   bool WantNLComments() const { return true; }

   const char* ObjDescription(int i)
   {  return csa->row_names[csa->obj_row[i]];
   }

   int ObjType(int i)
   {  int kind = mpl_get_row_kind(tran, csa->obj_row[i]+1);
      return kind == MPL_MIN ? NLW2_ObjSenseMinimize
         : NLW2_ObjSenseMaximize;
   }

   template <class ObjGradWriter>
   void FeedObjGradient(int i, ObjGradWriter& gw)
   {  assert(i < csa->m_obj);
      int *ind = csa->ind;
      double *val = csa->val;
      int len = mpl_get_mat_row(tran, csa->obj_row[i]+1, ind, val);
      if (len > 0)
      {  std::vector<std::pair<int, double>> tmp;
         tmp.reserve(len);
         for (int j = 1; j <= len; j++)
            if (ind[j] <= csa->n)
               tmp.emplace_back(csa->inv_perm[ind[j] - 1], val[j]);
         std::sort(tmp.begin(), tmp.end());
         auto sv = gw.MakeVectorWriter(tmp.size());
         for (auto const &p: tmp)
            sv.Write(p.first, p.second);
      }
   }

   template <class ExprWriter>
   void FeedObjExpression(int i, ExprWriter& ew)
   {  assert(i >= 0 && i < csa->m_obj);
      NL_EXPR *nl = mpl_get_row_nl(tran, csa->obj_row[i]+1);
      double c0 = mpl_get_row_c0(tran, csa->obj_row[i]+1);
      if (nl == NULL)
         ew.NPut(c0);
      else if (c0 == 0)
         convert_nl(ew, nl);
      else
      {  auto ew1 = ew.OPut2( mp::nl::ADD);
         ew1.NPut(c0);
         convert_nl(ew1, nl);
      }
   }

   template <class VarBoundsWriter>
   void FeedVarBounds(VarBoundsWriter& vbw)
   {  for (int i = 0; i < csa->n; i++)
      {  int j = csa->perm[i];
         int type, kind;
         double lb, ub;
         kind = mpl_get_col_kind(tran, j+1);
         type = mpl_get_col_bnds(tran, j+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (kind == MPL_BIN)
         {  if (type == MPL_FR || type == MPL_UP || lb < 0.0) lb = 0.0;
            if (type == MPL_FR || type == MPL_LO || ub > 1.0) ub = 1.0;
         }
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  if (fabs(lb) <= fabs(ub)) ub = lb;
            else lb = ub;
         }
         vbw.WriteLbUb(lb, ub);
      }
   }

   template <class IGWriter>
   void FeedInitialGuesses(IGWriter& igw)
   {  auto ig = igw.MakeVectorWriter(csa->n);
      for (int i = 0; i < csa->n; i++)
         ig.Write(i, tran->col[csa->perm[i]+1]->prim);
   }

   template <class ConBoundsWriter>
   void FeedConBounds(ConBoundsWriter& cbw)
   {  for (int i = 0; i < csa->m_con; i++)
      {  AlgConRange bnd;
         double lb, ub;
         int type;
         type = mpl_get_row_bnds(tran, csa->con_row[i]+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  if (fabs(lb) <= fabs(ub)) ub = lb;
            else lb = ub;
         }
         bnd.L = lb;
         bnd.U = ub;
         cbw.WriteAlgConRange(bnd);
      }
   }

   const char* ConDescription(int i) { return csa->row_names[csa->con_row[i]]; }

   template <class ConLinearExprWriter>
   void FeedLinearConExpr(int i, ConLinearExprWriter& clw)
   {  int *ind = csa->ind;
      double *val = csa->val;
      int len = mpl_get_mat_row(tran, csa->con_row[i]+1, ind, val);
      if (len > 0)
      {  std::vector<std::pair<int, double>> tmp;
         tmp.reserve(len);
         for (int j = 1; j <= len; j++)
            if (ind[j] <= csa->n)
               tmp.emplace_back(csa->inv_perm[ind[j] - 1], val[j]);
         // Some solvers crash if not in order
         std::sort(tmp.begin(), tmp.end());
         auto sv = clw.MakeVectorWriter(tmp.size());
         for (auto const &p: tmp)
            sv.Write(p.first, p.second);
      }
   }

   template <class ExprWriter>
   void FeedConExpression(int i, ExprWriter& ew)
   {  assert(i >= 0 && i < csa->m_con);
      NL_EXPR *nl = mpl_get_row_nl(tran, csa->con_row[i]+1);
      convert_nl(ew, nl);
   }

   template <class ColSizeWriter>
   void FeedColumnSizes(ColSizeWriter& csw)
   {  if (WantColumnSizes())
         for (int i = 0; i < csa->n - 1; i++)
            csw.Write(csa->col_size[csa->perm[i]]);
   }

   template <class RowObjNameWriter>
   void FeedRowAndObjNames(RowObjNameWriter& wrt)
   {  if (wrt)
      {  for (int i = 0; i < csa->m_con; i++)
            wrt << csa->row_names[csa->con_row[i]];
         for (int i = 0; i < csa->m_obj; i++)
            wrt << csa->row_names[csa->obj_row[i]];
      }
   }

   template <class ColNameWriter>
   void FeedColNames(ColNameWriter& wrt)
   {  if (wrt)
      {  for (int i = 0; i < csa->n; i++)
            wrt << csa->var_names[csa->perm[i]];
      }
   }

   template <class DefVarWriterFactory>
   void FeedDefinedVariables(int idx, DefVarWriterFactory& dvw)
   {  if (idx != 0) return;
      int n = tran->n;
      ELEMVAR **col = tran->col;
      for (int i = 0; i < csa->n_def; i++)
      {  int j = csa->defvar_order[i];
         std::vector<std::pair<int, double>> tmp;
         double c0 = 0;
         for (LINFORM *form = col[j+1]->def->lin; form != NULL; form = form->next)
         {  if (form->var == NULL)
            {  c0 = form->coef;
               continue;
            }
            tmp.emplace_back(csa->inv_perm[form->var->j-1], form->coef);
         }
         std::sort(tmp.begin(), tmp.end());
         auto dv = dvw.StartDefVar(j, tmp.size(), csa->var_names[j]);
         auto linw = dv.GetLinExprWriter();
         for (auto const &p: tmp)
            linw.Write(p.first, p.second);
         auto ew = dv.GetExprWriter();
         if (col[j+1]->def == NULL || col[j+1]->def->nl == NULL)
            ew.NPut(c0);
         else if (c0 == 0)
            convert_nl(ew, col[j+1]->def->nl);
         else
         {  auto ew1 = ew.OPut2(mp::nl::ADD);
            ew1.NPut(c0);
            convert_nl(ew1, col[j+1]->def->nl);
         }
      }
   }

private:

   const mp::nl::Opcode& get_opcode(int op)
   {  switch (op)
      {  case O_MUL: return mp::nl::MUL;
         case O_DIV: return mp::nl::DIV;
         case O_ADD: return mp::nl::ADD;
         case O_SUB: return mp::nl::SUB;
         case O_MINUS: return mp::nl::MINUS;
         case O_POWER: return mp::nl::POW;
         case O_SQRT: return mp::nl::SQRT;
         case O_ABS: return mp::nl::ABS;
         case O_COS: return mp::nl::COS;
         case O_SIN: return mp::nl::SIN;
         case O_TAN: return mp::nl::TAN;
         case O_ATAN: return mp::nl::ATAN;
         case O_ATAN2: return mp::nl::ATAN2;
         case O_EXP: return mp::nl::EXP;
         case O_LOG: return mp::nl::LOG;
         case O_LOG10: return mp::nl::LOG10;
         case O_SUM: return mp::nl::SUM;
         case O_MIN: return mp::nl::MIN;
         case O_MAX: return mp::nl::MAX;
         default:
            assert(op != op);
      }
   }

   template <class ExprWriter>
   void convert_nl(ExprWriter& ew, NL_EXPR *nl)
   {  if (nl == NULL)
      {  ew.NPut(0.0);
         return;
      }
      switch (nl->op)
      {  case O_MUL:
         case O_DIV:
         case O_ADD:
         case O_SUB:
         case O_POWER:
         case O_ATAN2:
            {  auto ew1 = ew.OPut2(get_opcode(nl->op));
               convert_nl(ew1, nl->u.arg.x);
               convert_nl(ew1, nl->u.arg.y);
            }
            break;
         case O_MINUS:
         case O_SQRT:
         case O_ABS:
         case O_COS:
         case O_SIN:
         case O_TAN:
         case O_ATAN:
         case O_EXP:
         case O_LOG:
         case O_LOG10:
            {  auto ew1 = ew.OPut1(get_opcode(nl->op));
               convert_nl(ew1, nl->u.arg.x);
            }
            break;
         case O_LINFORM:
            convert_linform(ew, nl->u.lin);
            break;
         case O_DEFVAR:
            ew.VPut(nl->u.defvar->j-1);
            break;
         case O_SUM:
         case O_MIN:
         case O_MAX:
            assert(nl->u.vec.n > 1);
            assert(!(nl->op == O_SUM && nl->u.vec.n == 2));
            {  auto ew1 = ew.OPutN(get_opcode(nl->op), nl->u.vec.n);
               for (int i = 0; i < nl->u.vec.n; i++)
                  convert_nl(ew1, nl->u.vec.data[i]);
            }
            break;
         default:
            assert(nl != nl);
      }
   }

   template <class ExprWriter>
   void convert_linform(ExprWriter& ew, LINFORM *lin)
   {  if (lin == NULL)
         ew.NPut(0.0);
      else if (lin->next == NULL)
         convert_linform_term(ew, lin);
#if 1
      else if (lin->next->next == NULL)
      {  auto ew1 = ew.OPut2(mp::nl::ADD);
         convert_linform_term(ew1, lin);
         convert_linform_term(ew1, lin->next);
      }
      else
      {  int ct = 0;
         LINFORM *t;
         for (t = lin; t; t = t->next)
            ct++;
         auto ew1 = ew.OPutN(mp::nl::SUM, ct);
         for (t = lin; t; t = t->next)
            convert_linform_term(ew1, t);
      }
#else
      else
      {  auto ew1 = ew.OPut2(mp::nl::ADD);
         convert_linform_term(ew1, lin);
         convert_linform(ew1, lin->next);
      }
#endif
   }

   template <class ExprWriter>
   void convert_linform_term(ExprWriter& ew, LINFORM *lin)
   {  if (lin->var == NULL)
      {  ew.NPut(lin->coef);
         return;
      }
      int j = lin->var->j-1;
      if (j < csa->n) j = csa->inv_perm[j];
      if (lin->coef == 1)
         ew.VPut(j);
      else if (lin->coef == -1)
      {  auto ew1 = ew.OPut1(mp::nl::MINUS);
         ew1.VPut(j);
      }
      else
      {  auto ew1 = ew.OPut2(mp::nl::MUL);
         ew1.NPut(lin->coef);
         ew1.VPut(j);
      }
   }
   mpl_tran *tran;
   struct csa *csa;
   mp::NLHeader hdr;
};


class Handler : public mp::SOLHandler
{
public:
   Handler(struct csa *c, const mp::NLHeader& h) : csa(c), hdr(h) {}

   mp::NLHeader Header() const { return hdr; }

   void OnSolveMessage(const char* s, int nbs)
   {  if (s!= NULL && strlen(s) > nbs)
         xprintf("%s\n", s + nbs);
   }

   int OnAMPLOptions(const AMPLOptions& ) { return 0; }

   template <class VecReader>
   void OnDualSolution(VecReader& rd)
   {  assert (rd.Size() == csa->m_con);
      csa->y = (double *)xcalloc(csa->m, sizeof(double));
      for (int i = 0; i < csa->m_con; i++)
         csa->y[csa->con_row[i]] = rd.ReadNext();
      for (int i = 0; i < csa->m_obj; i++)
         csa->y[csa->obj_row[i]] = 0;
   }

   template <class VecReader>
   void OnPrimalSolution(VecReader& rd)
   {  int n = rd.Size();
      assert (n == csa->n);
      csa->x = (double *)xcalloc(n, sizeof(double));
      for (int i = 0; rd.Size(); i++)
         csa->x[csa->perm[i]] = rd.ReadNext();
   }

   void OnObjno(int on) { csa->obj_no = on; }

   void OnSolveCode(int sr) { csa->solve_result = sr; }

   template <class SuffixReader>
   void OnIntSuffix(SuffixReader& sr)
   {  const auto& si = sr.SufInfo();
      int kind = si.Kind();
      const std::string& name = si.Name();
      assert(0 == (kind & 4));
      if (name == "sstatus")
      {  if (kind == 0)
         {  csa->var_status = (double *)xcalloc(csa->n, sizeof(double));
            for (int i = 0; i < csa->n; i++)
               csa->var_status[i] = 0;
            while (sr.Size())
            {  std::pair<int, int> val = sr.ReadNext();
               assert(val.first < csa->n);
               csa->var_status[csa->perm[val.first]] = val.second;
            }
         }
         else if (kind == 1)
         {  csa->row_status = (double *)xcalloc(csa->m, sizeof(double));
            for (int i = 0; i < csa->m; i++)
               csa->row_status[i] = 0;
            while (sr.Size())
            {  std::pair<int, int> val = sr.ReadNext();
               assert(val.first < csa->m_con);
               csa->row_status[csa->con_row[val.first]] = val.second;
            }
         }
      }
   }

private:
   struct csa *csa;
   mp::NLHeader hdr;
};

int mathprog_nlw2_run_solver(mpl_tran *tran, struct csa *csa,
   const char *solver, const char *options)
{     /* run solver and cleanup on error */
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_nlw2_run_solver: invalid call sequence\n");

      Feeder feeder(tran, csa);
      Handler handler(csa, feeder.Header());
      mp::NLSolver nlse;
      nlse.SetFileStub("");
      if (!nlse.Solve(feeder, handler, solver, options))
      {  xprintf("%s\n", nlse.GetErrorMessage());
         return 0;
      }
      return 1;
}

static double eval_linform(LINFORM *form, double *x)
{     double sum = 0;
      for (; form != NULL; form = form->next)
      {  if (form->var == NULL)
            sum += form->coef;
         else
            sum += form->coef * x[form->var->j-1];
      }
      return sum;
}

static double eval_nl_expr(NL_EXPR *nl, double *x)
{     if (nl == NULL)
         return 0;
      switch (nl->op)
      {  case O_MUL:
            return eval_nl_expr(nl->u.arg.x, x)
               * eval_nl_expr(nl->u.arg.y, x);
         case O_DIV:
            return eval_nl_expr(nl->u.arg.x, x)
               / eval_nl_expr(nl->u.arg.y, x);
         case O_ADD:
            return eval_nl_expr(nl->u.arg.x, x)
               + eval_nl_expr(nl->u.arg.y, x);
         case O_SUB:
            return eval_nl_expr(nl->u.arg.x, x)
               - eval_nl_expr(nl->u.arg.y, x);
         case O_MINUS:
            return -eval_nl_expr(nl->u.arg.x, x);
         case O_POWER:
            return pow(eval_nl_expr(nl->u.arg.x, x),
               eval_nl_expr(nl->u.arg.y, x));
         case O_SQRT:
            return sqrt(eval_nl_expr(nl->u.arg.x, x));
         case O_ABS:
            return fabs(eval_nl_expr(nl->u.arg.x, x));
         case O_COS:
            return cos(eval_nl_expr(nl->u.arg.x, x));
         case O_SIN:
            return sin(eval_nl_expr(nl->u.arg.x, x));
         case O_TAN:
            return tan(eval_nl_expr(nl->u.arg.x, x));
         case O_ATAN:
            return atan(eval_nl_expr(nl->u.arg.x, x));
         case O_ATAN2:
            return atan2(eval_nl_expr(nl->u.arg.x, x),
               eval_nl_expr(nl->u.arg.y, x));
         case O_EXP:
            return exp(eval_nl_expr(nl->u.arg.x, x));
         case O_LOG:
            return log(eval_nl_expr(nl->u.arg.x, x));
         case O_LOG10:
            return log10(eval_nl_expr(nl->u.arg.x, x));
         case O_LINFORM:
            return eval_linform(nl->u.lin, x);
         case O_DEFVAR:
            {  double sum;
               if (nl->u.defvar->def == NULL) return 0;
               sum = eval_linform(nl->u.defvar->def->lin, x);
               sum += eval_nl_expr(nl->u.defvar->def->nl, x);
               return sum;
            }
         case O_SUM:
            {  double sum = 0;
               for (int i = 0; i < nl->u.vec.n; i++)
                  sum += eval_nl_expr(nl->u.vec.data[i], x);
               return sum;
            }
         case O_MIN:
            {  double min = std::numeric_limits<double>::infinity();
               for (int i = 0; i < nl->u.vec.n; i++)
                  min = std::min(min, eval_nl_expr(nl->u.vec.data[i], x));
               return min;
            }
         case O_MAX:
            {  double max = -std::numeric_limits<double>::infinity();
               for (int i = 0; i < nl->u.vec.n; i++)
                  max = std::max(max, eval_nl_expr(nl->u.vec.data[i], x));
               return max;
            }
         default:
            assert(nl != nl);
      }
      assert(nl != nl);
      return 0;
}

int mathprog_nlw2_postsolve(mpl_tran *tran, struct csa *csa, int relax)
{     /* postsolve the model */
      int i, j, m, n, stat, ret, kind, len, t;
      double prim, dual;
      int *ind;
      double *val, *dval;

      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_nlw2_postsolve: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);

      ind = csa->ind;
      val = csa->val;
      dval = (double *)xcalloc(n, sizeof(double));
      for (j = 0; j < n; j++) dval[j] = 0;

      int mip = 0;
      if (!relax)
         for (j=0; j<n; j++)
            if (mpl_get_col_kind(tran, j+1) != MPL_NUM)
            {  mip = 1;
               break;
            }

      /* initialize var dual */
      if (!mip && csa->y != NULL && csa->obj_no >= 0)
      {  assert(csa->obj_no < csa->m_obj);
         len = mpl_get_mat_row(tran, csa->obj_row[csa->obj_no]+1, ind, val);
         for (t = 1; t <= len; t++)
            dval[ind[t]-1] = val[t];
      }

      for (i = 0; i < m; i++)
      {  prim = 0;
         stat = 0;
         dual = 0;
         len = mpl_get_mat_row(tran, i+1, ind, val);
         if (csa->x != NULL)
         {  prim = mpl_get_row_c0(tran, i+1);
            for (t = 1; t <= len; t++)
               prim += csa->x[ind[t]-1] * val[t];
            prim += eval_nl_expr(mpl_get_row_nl(tran, i+1), csa->x);
         }
         if (csa->y != NULL)
         {  dual = csa->y[i];
            if (!mip && csa->obj_no >= 0)
               for (t = 1; t <= len; t++)
                  dval[ind[t]-1] -= val[t] * csa->y[i];
         }
         stat = csa->row_status != NULL ? csa->row_status[i] : 0;
         if (stat == 2 || stat == 6) stat -= 2;
         if (stat == 3 || stat == 4) stat--;
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_row_soln(tran, i+1, stat, prim, dual);
      }

      for (j = 0; j < n; j++)
      {  prim = csa->x != NULL ? csa->x[j] : 0;
         stat = csa->var_status != NULL ? csa->var_status[j] : 0;
         dual = !mip ? dval[j] : 0;
         if (stat == 2 || stat == 6) stat = 0;
         if (stat == 3 || stat == 4) stat--;
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_col_soln(tran, j+1, stat, prim, dual);
      }

      ret = mpl_postsolve(tran);
      if (ret == 3) ret = csa->solve_result;
      else ret = -1;

      /* free working arrays */
      xfree(dval);
      return ret;
}

void mathprog_nlw2_cleanup(struct csa *csa)
{     int i, j;
      for (i = 0; i < csa->m; i++) xfree(csa->row_names[i]);
      xfree(csa->row_names);
      for (j = 0; j < csa->n + csa->n_def; j++) xfree(csa->var_names[j]);
      xfree(csa->var_names);
      xfree(csa->col_size);
      xfree(csa->prob_name);
      xfree(csa->obj_row);
      xfree(csa->con_row);

      if (csa->x)
         xfree(csa->x);
      if (csa->y)
         xfree(csa->y);
      if (csa->var_status)
         xfree(csa->var_status);
      if (csa->row_status)
         xfree(csa->row_status);
      xfree(csa->perm);
      xfree(csa->inv_perm);
      xfree(csa->ind);
      xfree(csa->val);
      if (csa->defvar_order)
         xfree(csa->defvar_order);

      delete csa;
}

#endif

/* eof */

