# TEST: category=exec_external
# TEST: solver=highs
# TEST: expect_status=0

var x >= 0;
maximize obj: x;
s.t. cap: x <= 2;

end;
