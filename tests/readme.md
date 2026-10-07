## Running tests

GoogleTest-based tests are wired through Meson.

```bash
meson setup builddir
meson compile -C builddir
meson test -C builddir
```

By default, external-solver execution tests are skipped. To enable them:

```bash
MATHPROG_ENABLE_EXTERNAL_TESTS=1 meson test -C builddir
```

Optionally override the external solver name declared in tests:

```bash
MATHPROG_ENABLE_EXTERNAL_TESTS=1 MATHPROG_EXTERNAL_SOLVER=<solver> meson test -C builddir
```
