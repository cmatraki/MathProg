# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=3.0
# TEST: obj_tol=1e-7

var x >= 0;
var y >= 0;
maximize obj: 3 * x + 2 * y;
s.t. cap: x + y <= 1;

end;
