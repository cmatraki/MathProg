This is a standalone extended version of GNU Mathprog, extracted from the GLPK package.
It is based on the fork by Domingo Alvarez Duarte (https://github.com/mingodad/GLPK) that
introduced several major user visible features including:
* let statements
* if/then/else and repeat statements
* break/continue statements
* sqlite3 driver for tables
* parsing of option statements
* named problems and multiple solves
* nested input files

This version adds:
* interfacing with solvers using .nl input (via nl-writer2)
* nonlinear and logical expressions
* defined variables (common subexpressions)
* fix/unfix statements
* actual option handling
* repeat until/while statements

## Bulding Mathprog

Mathprob is built using meson. For a quick start you can use the
following with any builddir.

```shell
meson setup builddir
meson compile -C builddir
```

Then you can run a model (say test.mod) using

```shell
builddir/mathprog [-m] test.mod [-d test.dat]
```

You can also run the test suite with

```shell
meson test -C builddir
```

See the file COPYING for the GNU General Public License.


