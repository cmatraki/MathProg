# TEST: category=exec_glpk
# TEST: expect_status=200

# LP HAS NO PRIMAL FEASIBLE SOLUTION

var x1;
var x2;

minimize obj: x1 + x2;

c1: x1 - x2 >= 5;
c2: -x1 + x2 >= 5;

solve;

end;
