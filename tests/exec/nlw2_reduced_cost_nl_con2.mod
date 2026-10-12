# TEST: category=exec_external
# TEST: solver=ipopt
# TEST: expect_status=0

var x >= 0;
var y >= 0;

minimize obj: 5 * x + y;

s.t. nlcon: x + y * y >= 4;

solve;

check abs(x) <= 1e-6 and abs(y - 2) <= 1e-6;
check abs(abs(nlcon.dual) - 0.25) <= 1e-6;
check abs(abs(x.dual) - 4.75) <= 1e-6 and abs(abs(y.dual)) <= 1e-6;

end;
