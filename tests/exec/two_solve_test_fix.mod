# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=0
# TEST: obj_tol=1e-7

var x >= 3;
var y >= -3;

fix y := 0;

minimize obj: x + y;

solve;

check x = 3 and y = 0;

unfix y;

fix x;

display x, y;

solve;

check x = 3 and y = -3;

display x, y;

end;
