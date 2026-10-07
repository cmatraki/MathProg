# TEST: category=syntax
# TEST: expect_parse=success

var x >= 0;
maximize obj: x;
s.t. cap: x <= 1;

end;
