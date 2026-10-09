# TEST: category=exec_glpk
# TEST: expect_status=200

# PROBLEM HAS NO PRIMAL FEASIBLE SOLUTION

var x integer >= 5;
var y;

minimize obj: x - 2 * y;

subj to c1: -x + 21 * y >=  2;
subj to c2: -3 * x + 2 * y <= 1;
subj to c3: 20 * x + y <= 20;

solve;

end;
