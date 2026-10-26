# TEST: category=syntax
# TEST: expect_parse=success

var x;
var y;

s.t. c1: x >= 3;
s.t. c2: y <= -3;

problem p1: x, y, c1, c2;

minimize obj: x - y;

end;
