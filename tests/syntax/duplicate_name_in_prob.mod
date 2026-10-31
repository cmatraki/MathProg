# TEST: category=syntax
# TEST: expect_parse=error

var x;

s.t. c0: x = 0;

problem p: c0, c0;

end
