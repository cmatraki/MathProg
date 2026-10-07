# TEST: category=exec_glpk
# TEST: expect_status=200

var x, integer, >= 0, <= 1;
maximize obj: x;
s.t. impossible: x >= 2;

end;
