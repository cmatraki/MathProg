# TEST: category=syntax
# TEST: expect_parse=success

param x default 0;
param y = 2 * x;

check y = 0;
check;

let x := 1;

check y = 2;

end;
