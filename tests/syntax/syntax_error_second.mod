# TEST: category=syntax
# TEST: expect_parse=error

set S;
param a{S} # missing semicolon (intentional syntax error)
var x{s in S} >= a[s];

end;
