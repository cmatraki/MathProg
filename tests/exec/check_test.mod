# TEST: category=exec_glpk
# TEST: expect_status=0

var x >= 3;

let x := 5;
 
check x = 5;              # check the assigned value

check;                     # run the first check

check x = 4;              # check future value

let x := 4;

minimize obj: x;

solve;                     # run the second check

check abs(x - 3) < 1e-6;   # check the result

end;                       # run the final check
