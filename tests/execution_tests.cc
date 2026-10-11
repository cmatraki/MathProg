#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>

#include "mathprog.h"
#include "test_case_loader.h"

#ifndef MATHPROG_TEST_DIR
#define MATHPROG_TEST_DIR "tests"
#endif

namespace {

int meta_int(const TestCase& tc, const char* key, int fallback) {
  const auto it = tc.meta.find(key);
  if (it == tc.meta.end())
    return fallback;
  return std::atoi(it->second.c_str());
}

double meta_double(const TestCase& tc, const char* key, double fallback) {
  const auto it = tc.meta.find(key);
  if (it == tc.meta.end())
    return fallback;
  return std::atof(it->second.c_str());
}

std::string meta_string(const TestCase& tc, const char* key,
                        const char* fallback = "") {
  const auto it = tc.meta.find(key);
  if (it == tc.meta.end())
    return fallback;
  return it->second;
}

extern "C" struct cb_dat {
  mpl_solver_options options;
  int count;
  int rc;
};

extern "C" int solve_callback(mpl_tran *tran, int sol_type, void *udata) {
    struct cb_dat *cbd = (struct cb_dat *)udata;
    cbd->count++;
    cbd->options.sol_type = sol_type;
    cbd->rc = mathprog_solver_run(tran, &cbd->options);
    return cbd->rc;
}

}  // namespace

class ExecutionCaseTest : public ::testing::TestWithParam<TestCase> {};

TEST_P(ExecutionCaseTest, RunModel) {
  const TestCase& tc = GetParam();
  const std::string category = meta_string(tc, "category");

  if (category == "exec_external") {
#ifndef WITH_NLW2
    GTEST_SKIP() << "External solver tests require WITH_NLW2";
#else
    const char* enabled = std::getenv("MATHPROG_ENABLE_EXTERNAL_TESTS");
    if (enabled == nullptr || std::string(enabled) != "1")
      GTEST_SKIP() << "Set MATHPROG_ENABLE_EXTERNAL_TESTS=1 to run external solver tests";
#endif
  }

  mpl_tran* tran = mathprog_alloc_wksp();
  ASSERT_NE(tran, nullptr);

  struct cb_dat cbd;
  std::memset(&cbd, 0, sizeof(cbd));

#ifdef WITH_NLW2
  if (category == "exec_external") {
    cbd.options.use_nlw2 = 1;
    std::string solver = meta_string(tc, "solver");
    const char* env_solver = std::getenv("MATHPROG_EXTERNAL_SOLVER");
    if (env_solver != nullptr && *env_solver != '\0')
      solver = env_solver;
    cbd.options.nlw2_solver = solver.empty() ? nullptr : solver.c_str();
  }
#endif

  mathprog_set_solve_callback(tran, solve_callback, &cbd);
  mathprog_init_rand(tran, 1);
  ASSERT_EQ(mathprog_read_model(tran, tc.path.c_str(), 0), 0) << tc.path;
  ASSERT_EQ(mathprog_generate(tran, nullptr), 0) << tc.path;
  const int rc = cbd.count > 0 ? cbd.rc : mathprog_solver_run(tran, &cbd.options);

  const int expected_status = meta_int(tc, "expect_status", 0);
  EXPECT_EQ(rc, expected_status) << tc.path;
  if (tc.meta.find("expect_obj") != tc.meta.end()) {
    const double expect_obj = meta_double(tc, "expect_obj", 0.0);
    const double tol = meta_double(tc, "obj_tol", 1e-7);
    EXPECT_NEAR(mathprog_eval_string(tran, "obj"), expect_obj, tol) << tc.path;
  }

  mathprog_free_wksp(tran);
}

INSTANTIATE_TEST_SUITE_P(
    ExecutionTests,
    ExecutionCaseTest,
    ::testing::ValuesIn(load_test_cases(MATHPROG_TEST_DIR, "exec_")),
    [](const ::testing::TestParamInfo<ExecutionCaseTest::ParamType>& info) {
      return test_case_name(info.param);
    });
