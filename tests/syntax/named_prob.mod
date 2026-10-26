# TEST: category=syntax
# TEST: expect_parse=success

var x;
var y;

minimize obj: x + y;

s.t. c1: x >= 3;
s.t. c2: y >= -3;

problem p1: obj, x, y, c1, c2;

end;
