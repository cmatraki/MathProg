#ifdef WITH_NLW2

extern "C" {
#include "mpl.h"

typedef struct NLW2_NLSolver_C mlp_nlw2_solver;
extern mlp_nlw2_solver *mathprog_nlw2_build_prob(mpl_tran *tran, int relax);
extern int mathprog_nlw2_postsolve(mpl_tran *tran, mlp_nlw2_solver *nsol, int relax);
extern int mathprog_nlw2_run_solver(mpl_tran *tran, mlp_nlw2_solver *nlse, const char *solver, const char *options);
}



#include "api/c/nl-solver-c.h"

typedef struct NLW2_NLSolver_C mlp_nlw2_solver;

mlp_nlw2_solver *mathprog_nlw2_build_prob(mpl_tran *tran, int relax)
{     /* build LP/MIP problem instance from the model */
      int m, n, i, j, t, kind, type, len, *ind;
      double lb, ub, *val;
      char *prob_name, *obj_name, *name;
      double *var_lb, *var_ub, *row_lb, *row_ub, *obj_c;
      int *var_type;
      char **var_names, **row_names;
      size_t *A_start;
      int *A_index;
      double *A_value;
      NLW2_ObjSense obj_sense;
      double obj_c0;
      int nnz;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mpl_build_prob: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      ind = (int *)xcalloc(1+n, sizeof(int));
      val = (double *)xcalloc(1+n, sizeof(double));
      var_lb = (double *)xcalloc(n, sizeof(double));
      var_ub = (double *)xcalloc(n, sizeof(double));
      var_type = (int *)xcalloc(n, sizeof(int));
      var_names = (char **)xcalloc(n, sizeof(char *));
      row_lb = (double *)xcalloc(m, sizeof(double));
      row_ub = (double *)xcalloc(m, sizeof(double));
      row_names = (char **)xcalloc(m, sizeof(char *));
      A_start = (size_t *)xcalloc(m, sizeof(size_t));
      nnz = 0;
      for (i = 0; i < m; i++)
         nnz += mpl_get_mat_row(tran, i+1, NULL, NULL);
      A_index = (int *)xcalloc(nnz, sizeof(int));
      A_value = (double *)xcalloc(nnz, sizeof(double));
      obj_c = NULL;

      for (j = 0; j < n; j++)
      {  /* set column name */
         name = mpl_get_col_name(tran, j+1);
         var_names[j] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(var_names[j], name);
         kind = mpl_get_col_kind(tran, j+1);
         var_type[j] = kind == MPL_NUM || relax ? 0 : 1;

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
         var_lb[j] = lb;
         var_ub[j] = ub;
      }
      nnz = 0;
      for (i = 0; i < m; i++)
      {  name = mpl_get_row_name(tran, i+1);
         row_names[i] = (char *)xcalloc(strlen(name)+1, sizeof(char));
         strcpy(row_names[i], name);
         type = mpl_get_row_bnds(tran, i+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (type == MPL_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         row_lb[i] = lb;
         row_ub[i] = ub;
         A_start[i] = nnz;
         len = mpl_get_mat_row(tran, i+1, ind, val);
         A_start[i] = nnz;
         for (t = 1; t <= len; t++) {
            A_index[nnz] = ind[t] - 1;
            A_value[nnz] = val[t];
            nnz++;
         }
      }
      /* build objective function (the first objective is used) */
      for (i = 0; i < m; i++)
      {  kind = mpl_get_row_kind(tran, i+1);
         if (kind == MPL_MIN || kind == MPL_MAX)
         {  /* set objective name */
            name = mpl_get_row_name(tran, i+1);
            obj_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
            strcpy(obj_name, name);
            /* prepare dense row */
            obj_c = (double *)xcalloc(n, sizeof(double));
            for (j = 0; j < n; j++) obj_c[j] = 0;
            /* set optimization direction */
            obj_sense = kind == MPL_MIN ? NLW2_ObjSenseMinimize : NLW2_ObjSenseMaximize;
            /* set constant term */
            obj_c0 = mpl_get_row_c0(tran, i+1);
            /* set objective coefficients */
            len = mpl_get_mat_row(tran, i+1, ind, val);
            for (t = 1; t <= len; t++)
               obj_c[ind[t]-1] = val[t];
            break;
         }
      }
      name = mpl_get_prob_name(tran);
      prob_name = (char *)xcalloc(strlen(name)+1, sizeof(char));
      strcpy(prob_name, name);
      NLW2_NLModel_C nlme = NLW2_MakeNLModel_C(prob_name);
      NLW2_SetCols_C(&nlme, n, var_lb, var_ub, var_type);
      NLW2_SetColNames_C(&nlme, (const char *const *)var_names);
      NLW2_SetRows_C(&nlme, m, row_lb, row_ub,
                     NLW2_MatrixFormatRowwise, nnz,
                     A_start, A_index, A_value);
      NLW2_SetRowNames_C(&nlme, (const char *const *)row_names);
      if(obj_c)
         NLW2_SetLinearObjective_C(&nlme, obj_sense, obj_c0, obj_c);
      NLW2_SetObjName_C(&nlme, obj_name);
      NLW2_NLOptionsBasic_C nlopts = NLW2_MakeNLOptionsBasic_C_Default();
      nlopts.n_text_mode_ = 1;
      nlopts.want_nl_comments_ = 1;
      NLW2_NLSolver_C *nlse = (NLW2_NLSolver_C *)xcalloc(1, sizeof(NLW2_NLSolver_C));
      *nlse = NLW2_MakeNLSolver_C(NULL);
      NLW2_SetNLOptions_C(nlse, nlopts);
      NLW2_SetFileStub_C(nlse, "");
      NLW2_LoadNLModel_C(nlse, &nlme);
      NLW2_DestroyNLModel_C(&nlme);

      /* free working arrays */
      for (i = 0; i < m; i++) xfree(row_names[i]);
      xfree(row_names);
      xfree(row_ub);
      xfree(row_lb);
      for (j = 0; j < n; j++) xfree(var_names[j]);
      xfree(var_names);
      xfree(var_type);
      xfree(var_ub);
      xfree(var_lb);
      if(obj_c)
      {  xfree(obj_c);
         xfree(obj_name);
      }
      xfree(A_index);
      xfree(A_value);
      xfree(A_start);
      xfree(prob_name);
      xfree(ind);
      xfree(val);
      return nlse;
}

int mathprog_nlw2_run_solver(mpl_tran *tran, mlp_nlw2_solver *nlse, const char *solver, const char *options)
{     /* run solver and cleanup on error */
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_nlw2_run_solver: invalid call sequence\n");
      if (!NLW2_RunSolver_C (nlse, solver, options))
      {  xprintf("%s\n", NLW2_GetErrorMessage_C(nlse));
         NLW2_DestroyNLSolver_C(nlse);
         xfree(nlse);
         return 0;
      }
      return 1;
}

int mathprog_nlw2_postsolve(mpl_tran *tran, mlp_nlw2_solver *nlse, int relax)
{     /* postsolve the model */
      int i, j, m, n, stat, ret, kind, len, t;
      double prim, dual;
            int *ind;
            double *val, *dval;
            const double *var_status = NULL, *row_status = NULL;
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_nlw2_postsolve: invalid call sequence\n");

      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);

      NLW2_NLSolution_C sol = NLW2_ReadSolution_C(nlse);

      for (t = 0; t < sol.nsuf_; t++)
      {
         if (!strcmp("sstatus", sol.suffixes_[t].name_)) {
            if(sol.suffixes_[t].kind_ == 0)
               var_status = sol.suffixes_[t].values_;
            else if(sol.suffixes_[t].kind_ == 1)
               row_status = sol.suffixes_[t].values_;
         }
      }

      ind = (int *)xcalloc(1+n, sizeof(int));
      val = (double *)xcalloc(1+n, sizeof(double));
      dval = (double *)xcalloc(n, sizeof(double));
      for (j = 0; j < n; j++) dval[j] = 0;

      int mip = 0;
      if (!relax)
         for(j=0; j<n; j++)
            if (mpl_get_col_kind(tran, j+1) != MPL_NUM)
            {  mip = 1;
               break;
            }

      /* initialize var dual (the first objective is used) */
      if (!mip && sol.y_)
         for (i = 0; i < m; i++)
         {  kind = mpl_get_row_kind(tran, i+1);
            if ((kind != MPL_MIN && kind != MPL_MAX))
            {
               len = mpl_get_mat_row(tran, i+1, ind, val);
               for (t = 1; t <= len; t++)
                  dval[ind[t]-1] = val[t];
               break;
            }
         }

      for (i = 0; i < m; i++)
      {  prim = 0;
         stat = 0;
         dual = 0;
         if (sol.x_)
         {  prim = mpl_get_row_c0(tran, i+1);
            len = mpl_get_mat_row(tran, i+1, ind, val);
            for (t = 1; t <= len; t++)
               prim += sol.x_[ind[t]-1] * val[t];
         }
         if (sol.y_)
         {  dual = sol.y_[i];
            if (!mip)
               for (t = 1; t <= len; t++)
                  dval[ind[t]-1] -= val[t] * sol.y_[i];
         }
         stat = row_status != NULL ? row_status[i] : 0;
         if (stat == 2 || stat == 6) stat -= 2;
         if (stat == 3 || stat == 4) stat--;
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_row_soln(tran, i+1, stat, prim, dual);
      }

      for (j = 0; j < n; j++)
      {  prim = sol.x_ ? sol.x_[j] : 0;
         stat = var_status != NULL ? var_status[j] : 0;
         dual = !mip ? dval[j] : 0;
         if (stat == 2 || stat == 6) stat = 0;
         if (stat == 3 || stat == 4) stat--;
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_col_soln(tran, j+1, stat, prim, dual);
      }

      ret = mpl_postsolve(tran);
      if (ret == 3) ret = sol.solve_result_;
      else ret = -1;

      xfree(ind);
      xfree(val);
      xfree(dval);
      NLW2_DestroyNLSolver_C(nlse);
      xfree(nlse);

      return ret;
}

#endif

/* eof */

