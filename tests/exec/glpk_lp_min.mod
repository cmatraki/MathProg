# TEST: category=exec_glpk
# TEST: expect_status=0
# TEST: expect_obj=1.0
# TEST: obj_tol=1e-7

var x >= 1;
minimize obj: x;

end;
