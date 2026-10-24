# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=0
# TEST: obj_tol=1e-7

param a, default 1;
param b = 2 * a;

var x >= b;
var y >= -b;

minimize obj: x + y;

solve;

check x = 2 and y = -2;

display x, y;

let a := 2;

solve;

check x = 4 and y = -4;

display x, y;

end;
