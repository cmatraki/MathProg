#ifdef WITH_NLW2

#include "mp/nl-solver.h"
#include "mp/nl-feeder.h"
#include "mp/sol-handler.h"
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
      char *prob_name, *obj_name;
      char **var_names, **row_names;
      double *col_size;
      NLW2_ObjSense obj_sense;
      double obj_c0;
      size_t nnz;
      int solve_result;
      double *x, *y, *var_status, *row_status;
      int m_range, m_eqn, n_bin, n_int, nnz_obj;
      int *perm, *inv_perm;
      int obj_no;
      int obj_row;
      int *ind;
      double *val;
};

struct csa *mathprog_nlw2_build_prob(mpl_tran *tran, int relax)
{     /* build LP/MIP problem instance from the model */
      struct csa *csa = new struct csa;
      int m, n, i, j, t, kind, type, len, *ind, nnz;
      double lb, ub, *val;
      char *name;
      int *int_idx, *bin_idx;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mpl_build_prob: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      csa->ind = ind = (int *)xcalloc(1+n, sizeof(int));
      int_idx = (int *)xcalloc(n, sizeof(int));
      bin_idx = (int *)xcalloc(n, sizeof(int));
      csa->val = val = (double *)xcalloc(1+n, sizeof(double));
      csa->m = m;
      csa->n = n;
      csa->var_names = (char **)xcalloc(n, sizeof(char *));
      csa->row_names = (char **)xcalloc(m, sizeof(char *));
      csa->col_size = (double *)xcalloc(n, sizeof(double));
      csa->obj_name = NULL;
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
      csa->obj_sense = NLW2_ObjSenseMinimize;
      csa->obj_row = -1;
      csa->obj_c0 = 0;

      t = 0;
      for (j = 0; j < n; j++)
      {  name = mpl_get_col_name(tran, j+1);
         csa->var_names[j] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->var_names[j], name);
         kind = mpl_get_col_kind(tran, j+1);
         csa->col_size[j] = 0;
         /* split variables per kind for permutation array */
         if (kind == MPL_NUM || relax)
            csa->perm[t++] = j;
         else if (kind == MPL_BIN)
            bin_idx[csa->n_bin++] = j;
         else
            int_idx[csa->n_int++] = j;
      }
      /* the varialbles should be ordered continuous, binary, integer */
      memcpy(csa->perm + t, bin_idx, csa->n_bin * sizeof(int));
      memcpy(csa->perm + t + csa->n_bin, int_idx, csa->n_int * sizeof(int));
      /* initialize inverse permutation */
      for (j = 0; j < n; j++)
         csa->inv_perm[csa->perm[j]] = j;

      nnz = 0;
      for (i = 0; i < m; i++)
      {  name = mpl_get_row_name(tran, i+1);
         csa->row_names[i] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->row_names[i], name);
         type = mpl_get_row_bnds(tran, i+1, NULL, NULL);
         if (type == MPL_FX)
            csa->m_eqn++;
         else if (type != MPL_FR)
            csa->m_range++;
         len = mpl_get_mat_row(tran, i+1, ind, NULL);
         nnz += len;
         for (t = 1; t <= len; t++) {
            csa->col_size[ind[t]-1]++;
         }
      }
      csa->nnz = nnz;
      /* build objective function (the first objective is used) */
      for (i = 0; i < m; i++)
      {  kind = mpl_get_row_kind(tran, i+1);
         if (kind == MPL_MIN || kind == MPL_MAX)
         {  /* set objective name */
            name = mpl_get_row_name(tran, i+1);
            csa->obj_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
            strcpy(csa->obj_name, name);
            /* set optimization direction */
            csa->obj_sense = kind == MPL_MIN ? NLW2_ObjSenseMinimize
               : NLW2_ObjSenseMaximize;
            /* set constant term */
            csa->obj_c0 = mpl_get_row_c0(tran, i+1);
            /* set objective coefficients */
            len = mpl_get_mat_row(tran, i+1, NULL, NULL);
            csa->nnz_obj = len;
            csa->obj_row = i;
            break;
         }
      }
      name = mpl_get_prob_name(tran);
      csa->prob_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
      strcpy(csa->prob_name, name);

      /* free working arrays */
      xfree(bin_idx);
      xfree(int_idx);

      return csa;
}

class Feeder : public mp::NLFeeder<struct csa*, void*>
{
public:
   Feeder(mpl_tran *t, struct csa *c) : tran(t), csa(c)
   {  hdr.num_vars = csa->n;
      hdr.num_algebraic_cons = csa->m;
      hdr.num_objs = csa->obj_row >= 0 ? 1 : 0;
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

      // the following are not (yet?) supported by MathProg
      hdr.num_logical_cons = 0;
      hdr.num_nl_cons = 0;
      hdr.num_nl_objs = 0;
      hdr.num_compl_conds = 0;
      hdr.num_nl_compl_conds = 0;
      hdr.num_compl_dbl_ineqs = 0;
      hdr.num_compl_vars_with_nz_lb = 0;
      hdr.num_nl_net_cons = 0;
      hdr.num_linear_net_cons = 0;
      hdr.num_nl_vars_in_cons = 0;
      hdr.num_nl_vars_in_objs = 0;
      hdr.num_nl_vars_in_both = 0;
      hdr.num_linear_net_vars = 0;
      hdr.num_funcs = 0;
      hdr.num_nl_integer_vars_in_both = 0;
      hdr.num_nl_integer_vars_in_cons = 0;
      hdr.num_nl_integer_vars_in_objs = 0;
      hdr.num_common_exprs_in_both = 0;
      hdr.num_common_exprs_in_cons = 0;
      hdr.num_common_exprs_in_objs = 0;
      hdr.num_common_exprs_in_single_cons = 0;
      hdr.num_common_exprs_in_single_objs = 0;
   }

   mp::NLHeader Header() { return hdr; }

   bool WantNLComments() const { return true; }

   const char* ObjDescription(int i)
   {  return csa->obj_name ? csa->obj_name : "";
   }

   int ObjType(int ) { return csa->obj_sense; }

   template <class ObjGradWriter>
   void FeedObjGradient(int i, ObjGradWriter& gw)
   {  if (csa->obj_row >= 0)
      {  int *ind = csa->ind;
         double *val = csa->val;
         int len = mpl_get_mat_row(tran, csa->obj_row+1, ind, val);
         if (len > 0)
         {  auto sv = gw.MakeVectorWriter(len);
            for (int j = 1; j <= len; j++)
               sv.Write(csa->inv_perm[ind[j] - 1], val[j]);
         }
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
            if (fabs(lb) <= fabs(ub)) ub = lb;
            else lb = ub;
         vbw.WriteLbUb(lb, ub);
      }
   }

   template <class ConBoundsWriter>
   void FeedConBounds(ConBoundsWriter& cbw)
   {  for (int i = 0; i < csa->m; i++)
      {  AlgConRange bnd;
         double lb, ub;
         int type;
         type = mpl_get_row_bnds(tran, i+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         bnd.L = lb;
         bnd.U = ub;
         cbw.WriteAlgConRange(bnd);
      }
   }

   const char* ConDescription(int i) { return csa->row_names[i]; }

   template <class ConLinearExprWriter>
   void FeedLinearConExpr(int i, ConLinearExprWriter& clw)
   {  int *ind = csa->ind;
      double *val = csa->val;
      int len = mpl_get_mat_row(tran, i+1, ind, val);
      if (len > 0)
      {  auto sv = clw.MakeVectorWriter(len);
         for (int j = 1; j <= len; j++)
            sv.Write(csa->inv_perm[ind[j] - 1], val[j]);
      }
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
      {  for (int i = 0; i < csa->m; i++)
            wrt << csa->row_names[i];
         if (csa->obj_name)
            wrt << csa->obj_name;
      }
   }

   template <class ColNameWriter>
   void FeedColNames(ColNameWriter& wrt)
   {  if (wrt)
      {  for (int i = 0; i < csa->n; i++)
            wrt << csa->var_names[csa->perm[i]];
      }
   }

private:
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
   {  int m = rd.Size();
      assert (m == csa->m);
      csa->y = (double *)xcalloc(m, sizeof(double));
      for (int i = 0; rd.Size(); i++)
         csa->y[i] = rd.ReadNext();
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
               assert(val.first < csa->m);
               csa->row_status[val.first] = val.second;
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
      /* initialize var dual (the first objective is used) */
      if (!mip && csa->y != NULL)
         if (csa->obj_row > -1)
         {  len = mpl_get_mat_row(tran, csa->obj_row+1, ind, val);
            for (t = 1; t <= len; t++)
               dval[ind[t]-1] = val[t];
         }

      for (i = 0; i < m; i++)
      {  prim = 0;
         stat = 0;
         dual = 0;
         if (csa->x != NULL)
         {  prim = mpl_get_row_c0(tran, i+1);
            len = mpl_get_mat_row(tran, i+1, ind, val);
            for (t = 1; t <= len; t++)
               prim += csa->x[ind[t]-1] * val[t];
         }
         if (csa->y != NULL)
         {  dual = csa->y[i];
            if (!mip)
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
      for (j = 0; j < csa->n; j++) xfree(csa->var_names[j]);
      xfree(csa->var_names);
      if(csa->obj_name)
         xfree(csa->obj_name);
      xfree(csa->col_size);
      xfree(csa->prob_name);

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

      delete csa;
}

#endif

/* eof */

