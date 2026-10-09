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

#include "mathprog.h"
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

int mathprog_set_msg_lev(mpl_tran *tran, int msg_lev)
{
    int rc = tran->msg_lev;
    tran->msg_lev = msg_lev;
    return rc;
}

double mathprog_eval_string(mpl_tran *tran, const char *expr)
{     /* evaluate numerical expression in given string */
      return mpl_eval_string(tran, expr);
}

void mathprog_get_stats(mpl_tran *tran, mpl_stats *stats)
{     /* get problem statistics */
      stats->rows = mpl_get_num_rows(tran);
      stats->log = mpl_get_num_log(tran);
      stats->cols = mpl_get_num_cols(tran);
      stats->defs = mpl_get_num_defs(tran);
}

static int mathprog_build_prob(mpl_tran *tran, glp_prob *prob)
{     /* build LP/MIP problem instance from the model */
      int m, n, i, j, t, kind, type, len, *ind;
      double lb, ub, *val;
      if (tran->phase != GLP_TRAN_PHASE_GENERATE)
         xerror("mathprog_build_prob: invalid call sequence\n");
      /* erase the problem object */
      glp_erase_prob(prob);
      /* set problem name */
      glp_set_prob_name(prob, mpl_get_prob_name(tran));
      /* build rows (constraints) */
      m = mpl_get_num_rows(tran);
      if (mpl_get_num_log(tran) > 0)
      {  xprintf("mathprog_build_prob: logical constraints not supporte"
         "d (%s)\n", mpl_get_row_name(tran, m + 1));
         tran->flag_p = 1;
         return 1;
      }
      if (m > 0)
         glp_add_rows(prob, m);
      for (i = 1; i <= m; i++)
      {  /* check for nonlinearities */
         if (mpl_get_row_nl(tran, i) != NULL)
         {  xprintf("mathprog_build_prob: nonlinear expressions not sup"
               "ported (%s)\n", mpl_get_row_name(tran, i));
            tran->flag_p = 1;
            return 1;
         }
         /* set row name */
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

static int mathprog_postsolve(mpl_tran *tran, glp_prob *prob, int sol)
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
      return ret;
}


#ifdef WITH_NLW2
void *mathprog_nlw2_build_prob(mpl_tran *tran, int relax);
int mathprog_nlw2_run_solver(mpl_tran *tran, void *csa,
   const char *solver, const char *options);
int mathprog_nlw2_postsolve(mpl_tran *tran, void *csa, int relax);
void mathprog_nlw2_cleanup(void *csa);
#endif

int mathprog_solver_run(mpl_tran *tran, const mpl_solver_options *opt)
{     /* call the solver */
      int ret = 0, psol_type = GLP_SOL, status = GLP_UNDEF;
      int sol_type = opt->sol_type;
      glp_prob *prob;
      if (opt && opt->use_nlw2)
#if !defined(WITH_NLW2)
         xprintf("No nlw2 support, use_nlw2 option ignored\n");
#else
      {  const char *solver, *solver_opts, *cfg_opt;
         int relax = 0;
         void *nlse;

         cfg_opt = mpl_get_option(tran, "relax_integrality");
         relax = cfg_opt != NULL && atoi(cfg_opt) > 0;
         relax |= opt->nomip;
         if (sol_type == 1 || sol_type == 3)
         {  /*A_PROBLEM_LP or A_PROBLEM_IP*/
            xprintf("Need to set solver options to select simplex or int"
               "erior point method\n");
            relax = 1;
         }
         solver = opt->nlw2_solver;
         if (solver == NULL)
            solver = mpl_get_option(tran, "solver");
         if (solver == NULL)
         {  xprintf("No solver given\n");
            return -1;
         }
         solver_opts = opt->nlw2_solver_opts;
         if (solver_opts == NULL)
         {  const char *suffix = "_options";
            char *name = malloc(strlen(solver) + strlen(suffix) + 1);
            strcpy(name, solver);
            strcat(name, suffix);
            solver_opts = mpl_get_option(tran, name);
            free(name);
         }
         if (solver_opts == NULL)
            solver_opts = "";

         /* build the problem instance from the model */
         nlse = mathprog_nlw2_build_prob (tran, relax);
         if (!mathprog_nlw2_run_solver(tran, nlse, solver, solver_opts))
            ret = -1;
         else
            ret = mathprog_nlw2_postsolve(tran, nlse, relax);
         mathprog_nlw2_cleanup(nlse);
         return ret;
      }
#endif
      /* build the problem instance from the model */
      prob = glp_create_prob();
      if (mathprog_build_prob(tran, prob))
      {  ret = 500;
         goto done;
      }

      if (opt && opt->new_name != NULL)
         glp_set_prob_name(prob, opt->new_name);

      glp_sort_matrix(prob);

      if (opt && opt->hide)
      {  int i, j;
         glp_set_obj_name(prob, NULL);
         glp_delete_index(prob);
         for (i = glp_get_num_rows(prob); i >= 1; i--)
            glp_set_row_name(prob, i, NULL);
         for (j = glp_get_num_cols(prob); j >= 1; j--)
            glp_set_col_name(prob, j, NULL);
      }

      if (sol_type == 0) /*A_PROBLEM_AUTO*/
      {  if (glp_get_num_int(prob) + glp_get_num_bin(prob) > 0)
            sol_type = 2; /*A_PROBLEM_MIP*/
         else
            sol_type = 1; /*A_PROBLEM_LP*/
      }

      if (sol_type == 2 && opt && opt->nomip)
         sol_type = 1;

      switch (sol_type)
      {  case 1: /*A_PROBLEM_LP*/
            {  glp_smcp smcp;
               glp_init_smcp(&smcp);
               smcp.presolve = GLP_ON;
               if (opt && opt->quiet)
                  smcp.msg_lev = GLP_MSG_OFF;
               ret = glp_simplex(prob, &smcp);
            }
            status = glp_get_status(prob);
            psol_type = GLP_SOL;
            break;
         case 2: /*A_PROBLEM_MIP*/
            {  glp_iocp iocp;
               glp_init_iocp(&iocp);
               iocp.presolve = GLP_ON;
               if (opt && opt->quiet)
                  iocp.msg_lev = GLP_MSG_OFF;
               ret = glp_intopt(prob, &iocp);
            }
            status = glp_mip_status(prob);
            psol_type = GLP_MIP;
            break;
         case 3: /*A_PROBLEM_IP*/
            {  glp_iptcp iptcp;
               glp_init_iptcp(&iptcp);
               if (opt && opt->quiet)
                  iptcp.msg_lev = GLP_MSG_OFF;
               ret = glp_interior(prob, &iptcp);
            }
            status = glp_ipt_status(prob);
            psol_type = GLP_IPT;
            break;
         default:
            xprintf("Unknown solution type %d\n", sol_type);
            ret = 500;
            goto done;
      }
      switch (status)
      {  case GLP_OPT:    /* solution is optimal */
            ret = 0;      /* solved */
            break;
         case GLP_NOFEAS: /* no feasible solution exists */
            ret = 200;    /* infeasible */
            break;
         case GLP_UNBND:  /* solution is unbounded */
            if (glp_get_prim_stat(prob) == GLP_FEAS)
               ret = 300; /* unbounded, feasible solution returned */
            else
               ret = 350; /* unbounded, no feasible solution returned */
            break;
         case GLP_FEAS:   /* solution is feasible */
            if (ret == GLP_ETMLIM || ret == GLP_EITLIM || ret == GLP_EMIPGAP)
               ret = 400; /* limit */
            else
               ret = 100; /* solved? solution candidate returned but error
                             likely */
            break;
         case GLP_INFEAS: /* solution is infeasible */
         case GLP_UNDEF:  /* solution is undefined */
            if (ret == GLP_ENOPFS)
               ret = 200; /* infeasible */
            else if (ret == GLP_ENODFS)
               ret = 450; /* limit, problem is either infeasible or unbounded */
            else if (ret == GLP_EITLIM || ret == GLP_ETMLIM || ret == GLP_EMIPGAP)
               ret = 470; /* limit, no solution returned */
            else
               ret = 500; /* failure */
            break;
      }
      if (mathprog_postsolve(tran, prob, psol_type))
      {  xprintf("Error on postsolving model\n");
         ret = 500;
      }
done: glp_delete_prob(prob);
      return ret;
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

/* eof */
