# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=-4.0
# TEST: obj_tol=1e-7

var x >= -5;
var y := 1;

minimize obj: x + y;

s.t. c1: x >= 3;
s.t. c2: y >= -3;

# y is not mentioned so should be fixed at 1.
problem p1: obj, x, c2;

solve;

display x, y, obj, x+y;

check x = -5 and y = 1;

end;
