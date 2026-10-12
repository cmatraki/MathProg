# TEST: category=exec_external
# TEST: solver=ipopt
# TEST: expect_status=0

/* Regression test for nlw2 reduced-cost computation with nonlinear objective.
 *
 * Expected KKT behavior for:
 *   minimize (x - 3)^2
 *   subject to x >= 5
 *
 * Optimum: x = 5, objective gradient = 2 * (x - 3) = 4.
 * Reduced cost for x should therefore have magnitude 4
 * (sign convention may depend on solver/status mapping).
 */

var x >= 5;

minimize obj: (x - 3) * (x - 3);

solve;

check abs(x - 5) <= 1e-6;
check abs(abs(x.dual) - 4) <= 1e-6;

printf "nlw2_reduced_cost_qp: x=%g dual=%g\n", x, x.dual;

end;
