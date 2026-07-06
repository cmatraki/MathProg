/* mpl.c (processing model in GNU MathProg language) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2008-2016 Andrew Makhorin, Department for Applied
*  Informatics, Moscow Aviation Institute, Moscow, Russia. All rights
*  reserved. E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "mpl.h"
#include "glpk.h"

mpl_tran *mathprog_alloc_wksp(void)
{     /* allocate the MathProg translator workspace */
      mpl_tran *tran;
      tran = mpl_initialize();
      return tran;
}

void mathprog_init_rand(mpl_tran *tran, int seed)
{     /* initialize pseudo-random number generator */
      if (tran->phase != GLP_TRAN_PHASE_INITIAL)
         xerror("mathprog_init_rand: invalid call sequence\n");
      rng_init_rand(tran->rand, seed);
      return;
}

int mathprog_read_model(mpl_tran *tran, const char *fname, int skip)
{     /* read and translate model section */
      int ret;
      if (tran->phase != GLP_TRAN_PHASE_INITIAL)
         xerror("mathprog_read_model: invalid call sequence\n");
      ret = mpl_read_model(tran, (char *)fname, skip);
      if (ret == 1 || ret == 2)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      else
         xassert(ret != ret);
      return ret;
}

int mathprog_read_data(mpl_tran *tran, const char *fname)
{     /* read and translate data section */
      int ret;
      if (!(tran->phase == GLP_TRAN_PHASE_MODEL || tran->phase == GLP_TRAN_PHASE_DATA))
         xerror("mathprog_read_data: invalid call sequence\n");
      ret = mpl_read_data(tran, (char *)fname);
      if (ret == 2)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      else
         xassert(ret != ret);
      return ret;
}

int mathprog_generate(mpl_tran *tran, const char *fname)
{     /* generate the model */
      int ret;
      if (!(tran->phase == GLP_TRAN_PHASE_MODEL || tran->phase == GLP_TRAN_PHASE_DATA))
         xerror("mathprog_generate: invalid call sequence\n");
      ret = mpl_generate(tran, (char *)fname);
      if (ret == 3)
         ret = 0;
      else if (ret == 4)
         ret = 1;
      return ret;
}

void mathprog_set_solve_callback(mpl_tran *tran, mpl_solve_callback solve_cb, void *udata)
{
    tran->solve_callback = solve_cb;
    tran->solve_callback_udata = udata;
}

int mathprog_set_genall(mpl_tran *tran, int bflag)
{
    int rc = tran->gen_all;
    tran->gen_all = bflag;
    return rc;
}

int mathprog_set_show_delta(mpl_tran *tran, int bflag)
{
    int rc = tran->show_delta;
    tran->show_delta = bflag;
    return rc;
}

int mathprog_set_add_missing_param_values(mpl_tran *tran, int bflag)
{
    int rc = tran->add_missing_param_values;
    tran->add_missing_param_values = bflag;
    return rc;
}

int mathprog_set_msg_lev(mpl_tran *tran, int msg_lev)
{
    int rc = tran->msg_lev;
    tran->msg_lev = msg_lev;
    return rc;
}

const char *mathprog_get_option(mpl_tran *tran, const char *name)
{     /* get option string */
      return mpl_get_option(tran, name);
}

int mathprog_build_prob(mpl_tran *tran, glp_prob *prob)
{     /* build LP/MIP problem instance from the model */
      int m, n, i, j, t, kind, type, len, *ind;
      double lb, ub, *val;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mathprog_build_prob: invalid call sequence\n");
      /* erase the problem object */
      glp_erase_prob(prob);
      /* fail if the problem is not supported */
      if (tran->is_nl)
      {  xprintf("mathprog_build_prob: nonlinear expressions not suppor"
            "ted\n");
         tran->flag_p = 1;
         return 1;
      }
      /* set problem name */
      glp_set_prob_name(prob, mpl_get_prob_name(tran));
      /* build rows (constraints) */
      m = mpl_get_num_rows(tran);
      if (m > 0)
         glp_add_rows(prob, m);
      for (i = 1; i <= m; i++)
      {  /* set row name */
         glp_set_row_name(prob, i, mpl_get_row_name(tran, i));
         /* set row bounds */
         type = mpl_get_row_bnds(tran, i, &lb, &ub);
         switch (type)
         {  case MPL_FR: type = GLP_FR; break;
            case MPL_LO: type = GLP_LO; break;
            case MPL_UP: type = GLP_UP; break;
            case MPL_DB: type = GLP_DB; break;
            case MPL_FX: type = GLP_FX; break;
            default: xassert(type != type);
         }
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  type = GLP_FX;
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         }
         glp_set_row_bnds(prob, i, type, lb, ub);
         /* warn about non-zero constant term */
         if (mpl_get_row_c0(tran, i) != 0.0)
            xprintf("mathprog_build_prob: row %s; constant term %.12g i"
               "gnored\n",
               mpl_get_row_name(tran, i), mpl_get_row_c0(tran, i));
      }
      /* build columns (variables) */
      n = mpl_get_num_cols(tran);
      if (n > 0)
         glp_add_cols(prob, n);
      for (j = 1; j <= n; j++)
      {  /* set column name */
         glp_set_col_name(prob, j, mpl_get_col_name(tran, j));
         /* set column kind */
         kind = mpl_get_col_kind(tran, j);
         switch (kind)
         {  case MPL_NUM:
               break;
            case MPL_INT:
            case MPL_BIN:
               glp_set_col_kind(prob, j, GLP_IV);
               break;
            default:
               xassert(kind != kind);
         }
         /* set column bounds */
         type = mpl_get_col_bnds(tran, j, &lb, &ub);
         switch (type)
         {  case MPL_FR: type = GLP_FR; break;
            case MPL_LO: type = GLP_LO; break;
            case MPL_UP: type = GLP_UP; break;
            case MPL_DB: type = GLP_DB; break;
            case MPL_FX: type = GLP_FX; break;
            default: xassert(type != type);
         }
         if (kind == MPL_BIN)
         {  if (type == GLP_FR || type == GLP_UP || lb < 0.0) lb = 0.0;
            if (type == GLP_FR || type == GLP_LO || ub > 1.0) ub = 1.0;
            type = GLP_DB;
         }
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
         {  type = GLP_FX;
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         }
         glp_set_col_bnds(prob, j, type, lb, ub);
      }
      /* load the constraint matrix */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      for (i = 1; i <= m; i++)
      {  len = mpl_get_mat_row(tran, i, ind, val);
         glp_set_mat_row(prob, i, len, ind, val);
      }
      /* build objective function (the first objective is used) */
      for (i = 1; i <= m; i++)
      {  kind = mpl_get_row_kind(tran, i);
         if (kind == MPL_MIN || kind == MPL_MAX)
         {  /* set objective name */
            glp_set_obj_name(prob, mpl_get_row_name(tran, i));
            /* set optimization direction */
            glp_set_obj_dir(prob, kind == MPL_MIN ? GLP_MIN : GLP_MAX);
            /* set constant term */
            glp_set_obj_coef(prob, 0, mpl_get_row_c0(tran, i));
            /* set objective coefficients */
            len = mpl_get_mat_row(tran, i, ind, val);
            for (t = 1; t <= len; t++)
               glp_set_obj_coef(prob, ind[t], val[t]);
            break;
         }
      }
      /* free working arrays */
      xfree(ind);
      xfree(val);
      return 0;
}

int mathprog_postsolve(mpl_tran *tran, glp_prob *prob, int sol)
{     /* postsolve the model */
      int i, j, m, n, stat, ret, has_shift;
      double prim, dual;
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("mathprog_postsolve: invalid call sequence\n");
      if (!(sol == GLP_SOL || sol == GLP_IPT || sol == GLP_MIP))
         xerror("mathprog_postsolve: sol = %d; invalid parameter\n",
            sol);
      m = mpl_get_num_rows(tran);
      n = mpl_get_num_cols(tran);
      if (!(m == glp_get_num_rows(prob) &&
            n == glp_get_num_cols(prob)))
         xerror("mathprog_postsolve: wrong problem object\n");
      if (!mpl_has_solve_stmt(tran))
      {  ret = 0;
         goto done;
      }
      has_shift = glp_get_obj_coef(prob, 0) != 0.0;
      for (i = 1; i <= m; i++)
      {  if (sol == GLP_SOL)
         {  stat = glp_get_row_stat(prob, i);
            if(i == 1 && has_shift)
                prim = glp_get_obj_val(prob);
            else
                prim = glp_get_row_prim(prob, i);
            dual = glp_get_row_dual(prob, i);
         }
         else if (sol == GLP_IPT)
         {  stat = 0;
            if(i == 1 && has_shift)
                prim = glp_ipt_obj_val(prob);
            else
                prim = glp_ipt_row_prim(prob, i);
            dual = glp_ipt_row_dual(prob, i);
         }
         else if (sol == GLP_MIP)
         {  stat = 0;
            if(i == 1 && has_shift)
                prim = glp_mip_obj_val(prob);
            else
                prim = glp_mip_row_val(prob, i);
            dual = 0.0;
         }
         else
            xassert(sol != sol);
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_row_soln(tran, i, stat, prim, dual);
      }
      for (j = 1; j <= n; j++)
      {  if (sol == GLP_SOL)
         {  stat = glp_get_col_stat(prob, j);
            prim = glp_get_col_prim(prob, j);
            dual = glp_get_col_dual(prob, j);
         }
         else if (sol == GLP_IPT)
         {  stat = 0;
            prim = glp_ipt_col_prim(prob, j);
            dual = glp_ipt_col_dual(prob, j);
         }
         else if (sol == GLP_MIP)
         {  stat = 0;
            prim = glp_mip_col_val(prob, j);
            dual = 0.0;
         }
         else
            xassert(sol != sol);
         if (fabs(prim) < 1e-9) prim = 0.0;
         if (fabs(dual) < 1e-9) dual = 0.0;
         mpl_put_col_soln(tran, j, stat, prim, dual);
      }
      ret = mpl_postsolve(tran);
      if (ret == 3)
         ret = 0;
      else if (ret == 4)
         ret = 1;
done: return ret;
}

int mathprog_waiting_solve(mpl_tran *tran)
{
    return tran->stmt && tran->stmt->type == A_SOLVE;
}

void mathprog_free_wksp(mpl_tran *tran)
{     /* free the MathProg translator workspace */
      mpl_terminate(tran);
      return;
}

#ifdef WITH_NLW2
#include "api/c/nl-solver-c.h"

NLW2_NLSolver_C *glp_nlw2_mpl_build_prob(mpl_tran *tran, int relax)
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
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      var_lb = xcalloc(n, sizeof(double));
      var_ub = xcalloc(n, sizeof(double));
      var_type = xcalloc(n, sizeof(int));
      var_names = xcalloc(n, sizeof(char *));
      row_lb = xcalloc(m, sizeof(double));
      row_ub = xcalloc(m, sizeof(double));
      row_names = xcalloc(m, sizeof(char *));
      A_start = xcalloc(m, sizeof(size_t));
      nnz = 0;
      for (i = 0; i < m; i++)
         nnz += mpl_get_mat_row(tran, i+1, NULL, NULL);
      A_index = xcalloc(nnz, sizeof(int));
      A_value = xcalloc(nnz, sizeof(double));
      obj_c = NULL;

      for (j = 0; j < n; j++)
      {  /* set column name */
         name = mpl_get_col_name(tran, j+1);
         var_names[j] = xcalloc(strlen(name)+1, sizeof(char));
         strcpy(var_names[j], name);
         kind = mpl_get_col_kind(tran, j+1);
         var_type[j] = kind == MPL_NUM || relax ? 0 : 1;

         /* set column bounds */
         type = mpl_get_col_bnds(tran, j+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (kind == MPL_BIN)
         {  if (type == GLP_FR || type == GLP_UP || lb < 0.0) lb = 0.0;
            if (type == GLP_FR || type == GLP_LO || ub > 1.0) ub = 1.0;
         }
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
            if (fabs(lb) <= fabs(ub)) ub = lb; else lb = ub;
         var_lb[j] = lb;
         var_ub[j] = ub;
      }
      nnz = 0;
      for (i = 0; i < m; i++)
      {  name = mpl_get_row_name(tran, i+1);
         row_names[i] = xcalloc(strlen(name)+1, sizeof(char));
         strcpy(row_names[i], name);
         type = mpl_get_row_bnds(tran, i+1, &lb, &ub);
         if (type == MPL_FR || type == MPL_UP) lb = -INFINITY;
         if (type == MPL_FR || type == MPL_LO) ub = INFINITY;
         if (type == GLP_DB && fabs(lb - ub) < 1e-9 * (1.0 + fabs(lb)))
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
            obj_name = xcalloc(strlen(name)+1, sizeof(char));
            strcpy(obj_name, name);
            /* prepare dense row */
            obj_c = xcalloc(n, sizeof(double));
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
      prob_name = xcalloc(strlen(name)+1, sizeof(char));
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
      NLW2_NLSolver_C *nlse = xcalloc(1, sizeof(NLW2_NLSolver_C));
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

void glp_nlw2_mpl_cleanup(mpl_tran *tran, NLW2_NLSolver_C *nlse);

int glp_nlw2_mpl_postsolve(mpl_tran *tran, NLW2_NLSolver_C *nlse, int relax)
{     /* postsolve the model */
      int i, j, m, n, stat, ret, kind, len, t;
      double prim, dual;
            int *ind;
            double *val, *dval;
            const double *var_status = NULL, *row_status = NULL;
      if (!(tran->phase == GLP_TRAN_PHASE_GENERATE && !tran->flag_p))
         xerror("glp_nlw2_mpl_postsolve: invalid call sequence\n");

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

      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(double));
      dval = xcalloc(n, sizeof(double));
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
      glp_nlw2_mpl_cleanup(tran, nlse);

      return ret;
}

void glp_nlw2_mpl_cleanup(mpl_tran *tran, NLW2_NLSolver_C *nlse)
{     NLW2_DestroyNLSolver_C(nlse);
      tran->flag_p = 1;   /* in case we skipped postsolve */
      xfree(nlse);
}

#endif

/* eof */
