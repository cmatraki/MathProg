# TEST: category=exec_glpk
# TEST: expect_status=300

# LP HAS UNBOUNDED PRIMAL SOLUTION

var x;
var y;

minimize obj: x + y;

subj to c1: x - 2 * y >= 0;
subj to c2: -x + 2 * y <= 0;

option glpk_options 'pre:solve=0';

solve;

end;
