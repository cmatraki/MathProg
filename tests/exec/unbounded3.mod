# TEST: category=exec_glpk
# TEST: expect_status=300

# PROBLEM HAS NO DUAL FEASIBLE SOLUTION

var x;

minimize obj: x;

option glpk_options 'pre:solve=0';

solve;

end;
