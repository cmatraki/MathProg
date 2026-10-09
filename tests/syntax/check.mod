# TEST: category=syntax
# TEST: expect_parse=success

param x;
var y;

let x := 1;
let y := 2;

check x = 1;
check y = 2;

end;
