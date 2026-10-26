# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=-3.0
# TEST: obj_tol=1e-7

var x;

s.t. c1: x >= 3;

# c1 is not part of p1 and x should be fixed at 0.
problem p1;

var y := 1;

minimize obj: x + y;

s.t. c2: y >= -3;

solve;

display x, y, p1;

check x = 0 and y = -3;

end;
