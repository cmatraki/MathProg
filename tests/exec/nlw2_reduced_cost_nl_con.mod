# TEST: category=exec_external
# TEST: solver=ipopt
# TEST: expect_status=0

/* nlw2_reduced_cost_nl_con.mod - regression test for nonlinear constraint
   gradient contributions to reduced costs.

   Model:
     minimize  x
     s.t.      x^2 >= 4     (active at x = 2)
     var x >= 0;

   At the optimal x* = 2 (active constraint x^2 = 4):
     KKT stationarity: 1 - y * 2*x = 0  =>  y = 1/(2*2) = 0.25
     The nonlinear constraint gradient dg/dx = 2*x = 4 must be accounted
     for in the reduced-cost computation; without it, the dual would be
     miscomputed. */

var x >= 0;

minimize obj: x;
s.t. nlcon: x * x >= 4;

solve;

check abs(x - 2) <= 1e-6;
check abs(abs(nlcon.dual) - 0.25) <= 1e-6;

printf "nlw2_reduced_cost_nl_con: x=%g nlcon.dual=%g\n", x, nlcon.dual;

end;
