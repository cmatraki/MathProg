/* glpk.h */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000-2018 Andrew Makhorin, Department for Applied
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

#ifndef MATHPROG_H
#define MATHPROG_H

#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct mpl_tran mpl_tran;
/* MathProg translator workspace */

mpl_tran *mathprog_alloc_wksp(void);
/* allocate the MathProg translator workspace */

void mathprog_init_rand(mpl_tran *tran, int seed);
/* initialize pseudo-random number generator */

int mathprog_read_model(mpl_tran *tran, const char *fname, int skip);
/* read and translate model section */

int mathprog_read_data(mpl_tran *tran, const char *fname);
/* read and translate data section */

int mathprog_generate(mpl_tran *tran, const char *fname);
/* generate the model */

typedef int (*mpl_solve_callback)(mpl_tran *mpl, int sol_type, void *udata);
void mathprog_set_solve_callback(mpl_tran *tran, mpl_solve_callback solve_cb, void *udata);
/* set the solve callback that will be called when a solve
   statement is executed, the user_data pointer is passed
   back */

int mathprog_set_genall(mpl_tran *tran, int bflag);
/* set the gen_all option, pre generate all entities in the model
   return the previous value */

int mathprog_set_show_delta(mpl_tran *tran, int bflag);
/* set the show_delta option, show time/memory delta usage
   return the previous value */

int mathprog_set_msg_lev(mpl_tran *tran, int msg_lev);
/* set the message level (verbosity)
   return the previous value */

const char *mathprog_get_option(mpl_tran *tran, const char *name);
/* get option string */

int mathprog_build_prob(mpl_tran *tran, glp_prob *prob);
/* build LP/MIP problem instance from the model */

int mathprog_postsolve(mpl_tran *tran, glp_prob *prob, int sol);
/* postsolve the model */

int mathprog_waiting_solve(mpl_tran *tran);
/* check if the model execution is waiting for solve */

void mathprog_free_wksp(mpl_tran *tran);
/* free the MathProg translator workspace */

#ifdef __cplusplus
}
#endif

#endif

/* eof */
