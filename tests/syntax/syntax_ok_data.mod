# TEST: category=syntax
# TEST: expect_parse=success

param a;
var x >= 0;
maximize obj: a * x;

data;
param a := 2;

end;
