# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=3.0
# TEST: obj_tol=1e-7

var x, integer, >= 0, <= 2;
var y, integer, >= 0, <= 2;
maximize obj: x + y;
s.t. cap: 2 * x + y <= 4;

end;
