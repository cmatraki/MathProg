#ifdef WITH_NLW2

#include "mp/nl-solver.h"

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
      double *var_lb, *var_ub, *row_lb, *row_ub, *obj_c;
      int *var_type;
      char **var_names, **row_names;
      size_t *A_start;
      int *A_index;
      double *A_value;
      NLW2_ObjSense obj_sense;
      double obj_c0;
      size_t nnz;
      int solve_result;
      double *x, *y, *var_status, *row_status;
};

struct csa *mathprog_nlw2_build_prob(mpl_tran *tran, int relax)
{     /* build LP/MIP problem instance from the model */
      struct csa *csa = (struct csa *)xcalloc(1, sizeof(struct csa));
      int m, n, i, j, t, kind, type, len, *ind, nnz;
      double lb, ub, *val;
      char *name;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mpl_build_prob: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      ind = (int *)xcalloc(1+n, sizeof(int));
      val = (double *)xcalloc(1+n, sizeof(double));
      csa->m = m;
      csa->n = n;
      csa->var_lb = (double *)xcalloc(n, sizeof(double));
      csa->var_ub = (double *)xcalloc(n, sizeof(double));
      csa->var_type = (int *)xcalloc(n, sizeof(int));
      csa->var_names = (char **)xcalloc(n, sizeof(char *));
      csa->row_lb = (double *)xcalloc(m, sizeof(double));
      csa->row_ub = (double *)xcalloc(m, sizeof(double));
      csa->row_names = (char **)xcalloc(m, sizeof(char *));
      csa->A_start = (size_t *)xcalloc(m, sizeof(size_t));
      nnz = 0;
      for (i = 0; i < m; i++)
         nnz += mpl_get_mat_row(tran, i+1, NULL, NULL);
      csa->A_index = (int *)xcalloc(nnz, sizeof(int));
      csa->A_value = (double *)xcalloc(nnz, sizeof(double));
      csa->obj_c = NULL;
      csa->solve_result = -1;
      csa->x = NULL;
      csa->y = NULL;
      csa->var_status = NULL;
      csa->row_status = NULL;

      for (j = 0; j < n; j++)
      {  /* set column name */
         name = mpl_get_col_name(tran, j+1);
         csa->var_names[j] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->var_names[j], name);
         kind = mpl_get_col_kind(tran, j+1);
         csa->var_type[j] = kind == MPL_NUM || relax ? 0 : 1;

         /* set column bounds */
         type = mpl_get_col_bnds(tran, j+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (kind == MPL_BIN)
         {  if (type == MPL_FR || type == MPL_UP || lb < 0.0) lb = 0.0;
            if (type == MPL_FR || type == MPL_LO || ub > 1.0) ub = 1.0;
         }
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         csa->var_lb[j] = lb;
         csa->var_ub[j] = ub;
      }
      nnz = 0;
      for (i = 0; i < m; i++)
      {  name = mpl_get_row_name(tran, i+1);
         csa->row_names[i] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(csa->row_names[i], name);
         type = mpl_get_row_bnds(tran, i+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         csa->row_lb[i] = lb;
         csa->row_ub[i] = ub;
         csa->A_start[i] = nnz;
         len = mpl_get_mat_row(tran, i+1, ind, val);
         csa->A_start[i] = nnz;
         for (t = 1; t <= len; t++) {
            csa->A_index[nnz] = ind[t] - 1;
            csa->A_value[nnz] = val[t];
            nnz++;
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
            /* prepare dense row */
            csa->obj_c = (double *)xcalloc(n, sizeof(double));
            for (j = 0; j < n; j++) csa->obj_c[j] = 0;
            /* set optimization direction */
            csa->obj_sense = kind == MPL_MIN ? NLW2_ObjSenseMinimize
               : NLW2_ObjSenseMaximize;
            /* set constant term */
            csa->obj_c0 = mpl_get_row_c0(tran, i+1);
            /* set objective coefficients */
            len = mpl_get_mat_row(tran, i+1, ind, val);
            for (t = 1; t <= len; t++)
               csa->obj_c[ind[t]-1] = val[t];
            break;
         }
      }
      name = mpl_get_prob_name(tran);
      csa->prob_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
      strcpy(csa->prob_name, name);

      /* free working arrays */
      xfree(ind);
      xfree(val);

      return csa;
}

int mathprog_nlw2_run_solver(mpl_tran *tran, struct csa *csa,
   const char *solver, const char *options)
{     /* run solver and cleanup on error */
      int m = csa->m;
      int n = csa->n;
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_nlw2_run_solver: invalid call sequence\n");

      mp::NLModel nlme(csa->prob_name);
      nlme.SetCols({csa->n, csa->var_lb, csa->var_ub, csa->var_type});
      nlme.SetColNames(csa->var_names);
      nlme.SetRows(csa->m, csa->row_lb, csa->row_ub,
         {csa->m, NLW2_MatrixFormatRowwise, csa->nnz, csa->A_start,
         csa->A_index, csa->A_value});
      nlme.SetRowNames(csa->row_names);
      if(csa->obj_c)
         nlme.SetLinearObjective(csa->obj_sense, csa->obj_c0, csa->obj_c);
      nlme.SetObjName(csa->obj_name);
      auto nlopts = NLW2_MakeNLOptionsBasic_C_Default();
      nlopts.n_text_mode_ = 1;
      nlopts.want_nl_comments_ = 1;
      mp::NLSolver nlse;
      nlse.SetNLOptions(nlopts);
      nlse.SetFileStub("");
      nlse.LoadModel(static_cast<const mp::NLModel&>(nlme));

      if (!nlse.Solve(solver, options))
      {  xprintf("%s\n", nlse.GetErrorMessage());
         return 0;
      }

      auto sol = nlse.ReadSolution();

      if (!sol.x_.empty())
      {   csa->x = (double *)xcalloc(n, sizeof(double));
          memcpy(csa->x, sol.x_.data(), n*sizeof(double));
      }

      if (!sol.y_.empty())
      {   csa->y = (double *)xcalloc(m, sizeof(double));
          memcpy(csa->y, sol.y_.data(), m*sizeof(double));
      }

      for (const auto &suf : sol.suffixes_)
      {  if (suf.name_ == "sstatus")
         {  if (suf.kind_ == 0)
            {   csa->var_status = (double *)xcalloc(n, sizeof(double));
                memcpy(csa->var_status, suf.values_.data(), n*sizeof(double));
            }
            else if (suf.kind_ == 1)
            {   csa->row_status = (double *)xcalloc(m, sizeof(double));
                memcpy(csa->row_status, suf.values_.data(), m*sizeof(double));
            }
         }
      }

      csa->solve_result = sol.solve_result_;

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

      ind = (int *)xcalloc(1+n, sizeof(int));
      val = (double *)xcalloc(1+n, sizeof(double));
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
         for (i = 0; i < m; i++)
         {  kind = mpl_get_row_kind(tran, i+1);
            if ((kind != MPL_MIN && kind != MPL_MAX))
            {  len = mpl_get_mat_row(tran, i+1, ind, val);
               for (t = 1; t <= len; t++)
                  dval[ind[t]-1] = val[t];
               break;
            }
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
      xfree(ind);
      xfree(val);
      xfree(dval);
      return ret;
}

void mathprog_nlw2_cleanup(struct csa *csa)
{     int i, j;
      for (i = 0; i < csa->m; i++) xfree(csa->row_names[i]);
      xfree(csa->row_names);
      xfree(csa->row_ub);
      xfree(csa->row_lb);
      for (j = 0; j < csa->n; j++) xfree(csa->var_names[j]);
      xfree(csa->var_names);
      xfree(csa->var_type);
      xfree(csa->var_ub);
      xfree(csa->var_lb);
      if(csa->obj_c)
      {  xfree(csa->obj_c);
         xfree(csa->obj_name);
      }
      xfree(csa->A_index);
      xfree(csa->A_value);
      xfree(csa->A_start);
      xfree(csa->prob_name);

      if (csa->x)
         xfree(csa->x);
      if (csa->y)
         xfree(csa->y);
      if (csa->var_status)
         xfree(csa->var_status);
      if (csa->row_status)
         xfree(csa->row_status);

      xfree(csa);
}

#endif

/* eof */

