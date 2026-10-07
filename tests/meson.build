gtest_dep = dependency('gtest', main : true, fallback: ['gtest', 'gtest_main_dep'])

test_dir = meson.project_source_root() / 'tests'

t_syntax = executable('mathprog_syntax_tests',
  files('syntax_tests.cc', 'test_case_loader.cc'),
  include_directories : include_directories('../src'),
  dependencies : [mathprog_core_dep, gtest_dep],
  cpp_args : ['-DMATHPROG_TEST_DIR="' + test_dir + '"'],
  install : false,
)

t_exec = executable('mathprog_execution_tests',
  files('execution_tests.cc', 'test_case_loader.cc'),
  include_directories : include_directories('../src'),
  dependencies : [mathprog_core_dep, gtest_dep],
  cpp_args : ['-DMATHPROG_TEST_DIR="' + test_dir + '"'],
  install : false,
)

test('mathprog-syntax-tests', t_syntax)
test('mathprog-execution-tests', t_exec)
