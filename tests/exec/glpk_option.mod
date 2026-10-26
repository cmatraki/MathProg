# TEST: category=exec_glpk
# TEST: expect_status=0

var y;
var z;

minimize o: 0;

s.t. c: y = 0;

option glpk_options 'pre:solve=8';

end;
