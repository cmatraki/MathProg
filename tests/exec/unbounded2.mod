# TEST: category=exec_glpk
# TEST: expect_status=450

# PROBLEM HAS NO DUAL FEASIBLE SOLUTION

var x;

minimize obj: x;

solve;

end;
