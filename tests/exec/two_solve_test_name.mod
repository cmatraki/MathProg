# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=0
# TEST: obj_tol=1e-7

var x;
var y := 0;

minimize obj: x + y;

s.t. c1: x >= 3;
s.t. c2: y >= -3;

problem p1: x, c1;

solve;

check x = 3 and y = 0;

unfix y;

fix x;

display x, y;

problem p2: y, c2;

solve;

check x = 3 and y = -3;

display x, y;

end;
