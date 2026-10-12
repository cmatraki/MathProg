# TEST: category=exec_external
# TEST: solver=ipopt
# TEST: expect_status=0

var x >= 0;
var y >= 0;

minimize obj: 0.5*x*x + 5*x + y;

s.t. con: x + y >= 4;

solve;

check abs(x) <= 1e-6 and abs(y - 4) <= 1e-6;
check abs(abs(con.dual) - 1) <= 1e-6;
check abs(abs(x.dual) - 4) <= 1e-6 and abs(abs(y.dual)) <= 1e-6;

end;
