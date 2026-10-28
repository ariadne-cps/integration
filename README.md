# Ariadne Integration

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/integration/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/integration/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/integration/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/integration/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/integration/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/integration)

Ariadne Integration is the standalone C++20 integration layer used by Ariadne. It provides validated and approximate ODE integrators, flow bounders and inclusion integration on top of Ariadne Function.

## Layout

Public headers live under `include/integration`, implementations under `src`, C++ tests under `tests`, benchmarks under `benchmarks`, Python bindings and tests under `python`, and standalone examples under `tutorials`.

The umbrella header is `include/ariadne-integration.hpp`.

## Dependencies

Integration depends directly on:

- [ariadne-cps/function](https://github.com/ariadne-cps/function), included as Git submodule `ss/fun`.
- [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration), used for shared build configuration.

Function supplies Geometry, Algebra, and the lower Ariadne layers transitively. The configuration submodule is pinned to the same commit used by Function.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/integration.git
cd integration
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target everything
ctest --test-dir build --output-on-failure
```

## Python bindings

When a compatible Python development environment is available, the standalone build produces `pyariadne` from the Function bindings (including their transitive dependencies) plus the Integration bindings.

## Contribution guidelines

If you would like to contribute to Ariadne Integration, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>

## License

Ariadne Integration is released under the GNU General Public License v3.0.
