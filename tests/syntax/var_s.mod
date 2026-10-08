# TEST: category=syntax
# TEST: expect_parse=success

var s >= 0, <= 1;

# Check that s.* is not confused by s.t. parsing.
display s.val, s.dual, s.lb, s.ub;

end;
