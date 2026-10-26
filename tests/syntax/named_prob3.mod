# TEST: category=syntax
# TEST: expect_parse=success

var x;

s.t. c1: x >= 3;

problem p1: x, c1;

var y;

minimize obj: x - y;

s.t. c2: y <= -3;

display p1;

end;
