# Ariadne Integration

[![Format Status](https://github.com/ariadne-cps/integration/actions/workflows/format.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/format.yml)
[![macOS Status](https://github.com/ariadne-cps/integration/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/integration/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/integration/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/integration/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/integration/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/integration/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/integration/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/integration)

Ariadne Integration is the standalone C++20 integration layer used by Ariadne. It provides validated and approximate ODE integrators, flow bounders and inclusion integration on top of Ariadne Function.
Although independently maintained and usable, this repository is also part of the [Ariadne repository structure](https://github.com/ariadne-cps/ariadne/blob/main/doc/structure.md).

## Layout

Public headers live under `include/integration`, implementations under `src`, C++ tests under `tests`, benchmarks under `benchmarks`, Python bindings and tests under `python`, and standalone C++ and Python examples under `tutorial`.

The umbrella header is `include/ariadne-integration.hpp`.

## Dependencies

Integration depends directly on:

- [ariadne-cps/function](https://github.com/ariadne-cps/function), included as Git submodule `ss/fun`.
- [ariadne-cps/build](https://github.com/ariadne-cps/build), used for shared build configuration.
- [ariadne-cps/testing](https://github.com/ariadne-cps/testing), pinned as the test-only submodule `ss/tst` for fixture-local assertions and diagnostics.

Function supplies Geometry, Algebra, and the lower Ariadne layers transitively. The configuration submodule is pinned to the same commit used by Function. Unit tests use Ariadne Testing fixtures, including variants for each integrator and logging-verbosity setting.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/integration.git
cd integration
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target everything
ctest --test-dir build --output-on-failure
```

## Common validated step API

Include `ariadne-integration.hpp` and use `ValidatedStepIntegratorInterface::step`
with a `ValidatedIntegrationRequest`. Every successful call returns a
`ValidatedIntegrationStep` with physical start/end times and independently
certified `reach()` and `evolve()` boxes.

A request with a duration and a state bound asks for an exact step. For an
adaptive autonomous equation, pass `suggest(maximum_step)` without a state bound
(optionally with an absolute initial time). The same solver chooses its bound
and may reduce the step for accuracy. The adapter caps a longer returned flow
at the requested maximum; `time_step()` on the result is the certified duration.
`step_size_is_maximum()` and `has_bounding_box()` distinguish the requests;
accessing a missing bound throws `std::logic_error`. Adaptive continuation keeps
the previous endpoint map and uncertainty coordinates. Explicit-time and
inclusion requests still require exact steps and explicit bounds. Python exposes
the adaptive constructors as `ValidatedIntegrationRequest.with_maximum_step(...)`.

| Adapter | Numerical field | Structured result |
|---|---|---|
| `FlowStepIntegratorAdapter` | Single-valued, autonomous or explicit time with fixed parameters | Reach patches and endpoint map |
| `MultiflowStepIntegratorAdapter` | Box multifunction; timed fields require a timed backend | Boxes |
| `StructuredInclusionStepIntegratorAdapter` | Multifunction with the chosen backend's structural capability | Reach patches and endpoint map |

Construct the chosen backend and adapter explicitly. For structured inclusions,
configure exactly one reach approximation and one evolve approximation through
`SelectedParametricMultifunctionIntegrator`; they are never ranked or replaced
by an implicit fallback. `ZeroApproximation` and `ConstantApproximation` support
reach; the five approximation families remain available for endpoint evolution.

A continuation request takes the previous step and preserves its original state
coordinates, fixed parameters and accumulated input coordinates. Its time layout
continues the previous request; changing the field's argument layout is outside
this continuation contract. `reach_sections(absolute_time)` returns validated
state patches with the time coordinate removed. Use `has_reach_models()` and
`has_evolve_model()` before accessing structured output. A box-only backend
rejects a structured continuation instead of silently discarding correlations.
Reach and evolve can use different input coordinates; their domains need not be
identical. Declared fixed-parameter domains must match both maps.

Single-valued function handles are stored by value in requests, including
implicit conversions from effective functions. Polymorphic multifunction fields
are borrowed from named C++ objects, which must outlive the requests. Adapters
borrow their backends, which must outlive the adapters. Borrowing a temporary is
rejected at compilation. Python wrappers retain these dependencies automatically;
returned boxes and patches are independent values.

Unsupported capabilities or inconsistent domains raise `std::invalid_argument`
(`ValueError` in Python). A failed numerical certificate raises a flow exception
(`RuntimeError` in Python); it does not select another solver. The C++ and Python
[tutorials](tutorial) demonstrate two-step equation and inclusion integration,
structured output, and explicit failure handling.

## Python bindings

When a compatible Python development environment is available, import `pyariadne_integration`. It combines the Function bindings (including their transitive dependencies) and Integration bindings, including the common validated step API and the explicit first-order and structured inclusion adapters.

## Contribution guidelines

If you would like to contribute to Ariadne Integration, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>
