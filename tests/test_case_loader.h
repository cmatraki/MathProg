#ifndef TEST_CASE_LOADER_H
#define TEST_CASE_LOADER_H

#include <map>
#include <string>
#include <vector>

struct TestCase {
  std::string path;
  std::map<std::string, std::string> meta;
};

std::vector<TestCase> load_test_cases(const std::string& root,
                                      const std::string& category_prefix);
std::string test_case_name(const TestCase& tc);

#endif
