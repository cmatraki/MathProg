#include <gtest/gtest.h>

#include <string>
#include <vector>

extern "C" {
#include <glpk.h>
#include "mathprog.h"
}

#include "test_case_loader.h"

#ifndef MATHPROG_TEST_DIR
#define MATHPROG_TEST_DIR "tests"
#endif

class SyntaxCaseTest : public ::testing::TestWithParam<TestCase> {};

TEST_P(SyntaxCaseTest, ParseAndGenerate) {
  const TestCase& tc = GetParam();
  const auto it = tc.meta.find("expect_parse");
  const bool expect_success = (it == tc.meta.end() || it->second == "success");

  mpl_tran* tran = mathprog_alloc_wksp();
  ASSERT_NE(tran, nullptr);
  mathprog_init_rand(tran, 1);

  int ret = mathprog_read_model(tran, tc.path.c_str(), 0);
  if (ret == 0)
    ret = mathprog_generate(tran, nullptr);

  if (expect_success)
    EXPECT_EQ(ret, 0) << tc.path;
  else
    EXPECT_NE(ret, 0) << tc.path;

  mathprog_free_wksp(tran);
}

INSTANTIATE_TEST_SUITE_P(
    SyntaxTests,
    SyntaxCaseTest,
    ::testing::ValuesIn(load_test_cases(MATHPROG_TEST_DIR, "syntax")),
    [](const ::testing::TestParamInfo<SyntaxCaseTest::ParamType>& info) {
      return test_case_name(info.param);
    });
