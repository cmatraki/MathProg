project('MathProg', ['c', 'cpp'],
  version : '0.1.0',
  default_options : ['c_std=c99', 'cpp_std=c++17'])

cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

with_nlw2 = get_option('with_nlw2')
compile_defs = ['-DWITH_SPLAYTREE']
if with_nlw2
  compile_defs += ['-DWITH_NLW2']
endif

add_project_arguments(compile_defs, language : ['c', 'cpp'])

inc = include_directories('src', 'src/env', 'src/misc', 'src/mpl', 'nl-writer2/include')

glpk_dep = dependency('glpk', required : false)
if not glpk_dep.found()
  glpk_dep = cc.find_library('glpk', required : true)
endif
m_dep = cc.find_library('m', required : false)
z_dep = dependency('zlib', required : false)
if not z_dep.found()
  z_dep = cc.find_library('z', required : true)
endif
dl_dep = cc.find_library('dl', required : false)

deps = [glpk_dep, z_dep]
if m_dep.found()
  deps += m_dep
endif
if dl_dep.found()
  deps += dl_dep
endif

mathprog_core_sources = files(
  'src/api/mpl.c',
  'src/env/alloc.c',
  'src/env/dlsup.c',
  'src/env/env.c',
  'src/env/error.c',
  'src/env/stdc.c',
  'src/env/stream.c',
  'src/env/tls.c',
  'src/misc/avl.c',
  'src/misc/dmp.c',
  'src/misc/jd.c',
  'src/misc/splaytree.c',
  'src/misc/str2int.c',
  'src/misc/str2num.c',
  'src/misc/strspx.c',
  'src/misc/strtrim.c',
  'src/misc/rng.c',
  'src/misc/rng1.c',
  'src/mpl/mpl1.c',
  'src/mpl/mpl2.c',
  'src/mpl/mpl3.c',
  'src/mpl/mpl4.c',
  'src/mpl/mpl5.c',
  'src/mpl/mpl6.c',
  'src/mpl/mplsql.c',
)

if with_nlw2
  mathprog_core_sources += files(
    'src/mpl/mplnlw2.cc',
    'nl-writer2/src/nl-solver.cc',
    'nl-writer2/src/nl-utils.cc',
    'nl-writer2/src/nl-writer2.cc',
  )
endif

mathprog_core = static_library('mathprog_core',
  mathprog_core_sources,
  include_directories : inc,
  dependencies : deps,
  install : false,
)

mathprog_core_dep = declare_dependency(
  link_with : mathprog_core,
  include_directories : inc,
  dependencies : deps,
)

executable('mathprog',
  'src/mathprog.c',
  include_directories : inc,
  dependencies : [mathprog_core_dep],
  install : false,
)

if get_option('tests')
  subdir('tests')
endif
