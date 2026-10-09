# TEST: category=exec_glpk
# TEST: expect_status=450

# PROBLEM HAS NO DUAL FEASIBLE SOLUTION

var x >= 5;
var y;
var z;

minimize obj: x - 2 * y + 5 * z;

subj to c1: -x + 21 * y >=  2;
subj to c2: -3 * x + 2 * y <=  1;
subj to c3: 20 * x + y <= 20;

solve;

end;
