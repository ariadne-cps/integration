#include "pybind11.hpp"

void integration_submodule(pybind11::module& module);

void integration_module(pybind11::module& module) {
    integration_submodule(module);
}
