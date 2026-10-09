# TEST: category=syntax
# TEST: expect_parse=success

param x;
let x := 8;
repeat until x == 0 {
    check x > 0;
    check;
    let x := x-1;
}
check x = 0;

end;

