#include "test_case_loader.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>

namespace {

std::string trim(std::string s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
    s.erase(s.begin());
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
    s.pop_back();
  return s;
}

bool starts_with(const std::string& s, const std::string& prefix) {
  return s.size() >= prefix.size() &&
         std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::map<std::string, std::string> parse_meta(const std::string& path) {
  static const std::string kPrefix = "# TEST:";
  std::ifstream in(path);
  std::map<std::string, std::string> meta;
  std::string line;
  while (std::getline(in, line)) {
    if (!starts_with(line, kPrefix))
      continue;
    std::string body = trim(line.substr(kPrefix.size()));
    std::size_t eq = body.find('=');
    if (eq == std::string::npos)
      continue;
    std::string key = trim(body.substr(0, eq));
    std::string value = trim(body.substr(eq + 1));
    if (!key.empty())
      meta[key] = value;
  }
  return meta;
}

}  // namespace

std::vector<TestCase> load_test_cases(const std::string& root,
                                      const std::string& category_prefix) {
  std::vector<TestCase> out;
  namespace fs = std::filesystem;

  for (const auto& entry : fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file() || entry.path().extension() != ".mod")
      continue;

    TestCase tc;
    tc.path = entry.path().string();
    tc.meta = parse_meta(tc.path);

    auto it = tc.meta.find("category");
    if (it == tc.meta.end())
      continue;
    if (!starts_with(it->second, category_prefix))
      continue;

    out.push_back(std::move(tc));
  }

  std::sort(out.begin(), out.end(),
            [](const TestCase& a, const TestCase& b) { return a.path < b.path; });
  return out;
}

std::string test_case_name(const TestCase& tc) {
  std::filesystem::path p(tc.path);
  std::string stem = p.stem().string();
  for (char& ch : stem) {
    if (!std::isalnum(static_cast<unsigned char>(ch)))
      ch = '_';
  }
  return stem;
}
