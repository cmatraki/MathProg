# TEST: category=exec_glpk
# TEST: expect_status=470

var x >= 0;
var y >= 0;

maximize obj: 3 * x + 2 * y;

s.t. cap: x + y <= 1;

option glpk_options 'pre:solve=1 lim:iter=1';

solve;

end;
