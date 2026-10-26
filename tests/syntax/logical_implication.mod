# TEST: category=syntax
# TEST: expect_parse=error

param x;

let x := 0;

if x ==> 0 then let x:= 1;

end;
