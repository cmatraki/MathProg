# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=0.0
# TEST: obj_tol=1e-7

var x;

s.t. c1: x >= 3;

problem p1: x, c1;

var y := 1;

minimize obj: x + y;

s.t. c2: y >= -3;

solve;

display x, y, obj, x+y;

check x = 3 and y = -3;

end;
