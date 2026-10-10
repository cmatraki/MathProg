# TEST: category=exec_glpk
# TEST: expect_status=-1

var x >= 0;

minimize obj: x;

option glpk_options 'foo=bar';

solve;

end;
