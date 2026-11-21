/***************************************************************************
 *            validated_step_submodule.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/*
 *  This file is part of Ariadne Integration.
 *
 *  Ariadne Integration is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Integration is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file python/bindings/validated_step_submodule.cpp
 *  @brief Python access to the common validated numerical step contract.
 */

#include "function/taylor_function.hpp"
#include "integration/validated_step_integrator.hpp"
#include "pybind11.hpp"

using namespace Ariadne;

namespace {

std::vector<ValidatedVectorMultivariateFunctionPatch> python_patches(
    const List<ValidatedVectorMultivariateFunctionPatch>& patches) {
    return {patches.begin(), patches.end()};
}

} // namespace

Void export_validated_steps(pybind11::module& module) {
    namespace py = pybind11;

    py::class_<ValidatedIntegrationStep>(module, "ValidatedIntegrationStep")
        .def("time_step", &ValidatedIntegrationStep::time_step, py::return_value_policy::copy)
        .def("initial_time", &ValidatedIntegrationStep::initial_time, py::return_value_policy::copy)
        .def("final_time", &ValidatedIntegrationStep::final_time)
        .def("has_explicit_time", &ValidatedIntegrationStep::has_explicit_time)
        .def("constant_parameter_domain", &ValidatedIntegrationStep::constant_parameter_domain, py::return_value_policy::copy)
        .def("reach", &ValidatedIntegrationStep::reach, py::return_value_policy::copy)
        .def("evolve", &ValidatedIntegrationStep::evolve, py::return_value_policy::copy)
        .def("has_reach_models", &ValidatedIntegrationStep::has_reach_models)
        .def("reach_models", [](const ValidatedIntegrationStep& step) { return python_patches(step.reach_models()); })
        .def("reach_time_index", &ValidatedIntegrationStep::reach_time_index)
        .def("reach_time_offset", &ValidatedIntegrationStep::reach_time_offset, py::return_value_policy::copy)
        .def("reach_sections", [](const ValidatedIntegrationStep& step, const StepSizeType& time) {
            return python_patches(step.reach_sections(time));
        })
        .def("has_evolve_model", &ValidatedIntegrationStep::has_evolve_model)
        .def("evolve_model", &ValidatedIntegrationStep::evolve_model, py::return_value_policy::copy);

    py::class_<ValidatedBoxMultifunctionInterface>(module, "ValidatedBoxMultifunctionInterface")
        .def("argument_size", &ValidatedBoxMultifunctionInterface::argument_size)
        .def("result_size", &ValidatedBoxMultifunctionInterface::result_size)
        .def("domain", &ValidatedBoxMultifunctionInterface::domain)
        .def("evaluate", &ValidatedBoxMultifunctionInterface::evaluate);
    py::class_<ValidatedTimedBoxMultifunctionInterface>(module, "ValidatedTimedBoxMultifunctionInterface")
        .def("argument_size", &ValidatedTimedBoxMultifunctionInterface::argument_size)
        .def("result_size", &ValidatedTimedBoxMultifunctionInterface::result_size)
        .def("domain", &ValidatedTimedBoxMultifunctionInterface::domain)
        .def("time_domain", &ValidatedTimedBoxMultifunctionInterface::time_domain)
        .def("parameter_domain", &ValidatedTimedBoxMultifunctionInterface::parameter_domain)
        .def("evaluate", &ValidatedTimedBoxMultifunctionInterface::evaluate);
    py::class_<ValidatedInputFieldBoxAdapterDP, ValidatedBoxMultifunctionInterface>(module, "ValidatedInputFieldBoxAdapterDP")
        .def(py::init<const EffectiveVectorMultivariateFunction&, const BoxDomainType&, const BoxDomainType&>(),
             py::arg("field"), py::arg("inputs"), py::arg("states"));
    py::class_<ValidatedTimedInputFieldBoxAdapterDP, ValidatedTimedBoxMultifunctionInterface>(module, "ValidatedTimedInputFieldBoxAdapterDP")
        .def(py::init<const EffectiveVectorMultivariateFunction&, const BoxDomainType&, const IntervalDomainType&, const BoxDomainType&, const BoxDomainType&>(),
             py::arg("field"), py::arg("states"), py::arg("time"), py::arg("parameters"), py::arg("inputs"));

    py::class_<ValidatedIntegrationRequest> request(module, "ValidatedIntegrationRequest");
    request
        .def(py::init<const ValidatedVectorMultivariateFunction&, const ExactBoxType&, const StepSizeType&, const UpperBoxType&>())
        .def(py::init<const ValidatedVectorMultivariateFunction&, const ExactBoxType&, const StepSizeType&, const StepSizeType&, const ExactBoxType&, const UpperBoxType&>())
        .def(py::init<const ValidatedVectorMultivariateFunction&, const ValidatedIntegrationStep&, const StepSizeType&, const UpperBoxType&>())
        // Keep the Python field alive as long as this borrowing request.
        .def(py::init<const ValidatedBoxMultifunctionInterface&, const ExactBoxType&, const StepSizeType&, const UpperBoxType&>(), py::keep_alive<1, 2>())
        .def(py::init<const ValidatedBoxMultifunctionInterface&, const ValidatedIntegrationStep&, const StepSizeType&, const UpperBoxType&>(), py::keep_alive<1, 2>())
        .def(py::init<const ValidatedTimedBoxMultifunctionInterface&, const ExactBoxType&, const StepSizeType&, const StepSizeType&, const ExactBoxType&, const UpperBoxType&>(), py::keep_alive<1, 2>())
        .def(py::init<const ValidatedTimedBoxMultifunctionInterface&, const ValidatedIntegrationStep&, const StepSizeType&, const UpperBoxType&>(), py::keep_alive<1, 2>())
        .def("initial", &ValidatedIntegrationRequest::initial, py::return_value_policy::copy)
        .def("time_step", &ValidatedIntegrationRequest::time_step, py::return_value_policy::copy)
        .def("initial_time", &ValidatedIntegrationRequest::initial_time, py::return_value_policy::copy)
        .def("has_explicit_time", &ValidatedIntegrationRequest::has_explicit_time)
        .def("constant_parameter_domain", &ValidatedIntegrationRequest::constant_parameter_domain, py::return_value_policy::copy)
        .def("bounding_box", &ValidatedIntegrationRequest::bounding_box, py::return_value_policy::copy)
        .def("has_initial_mapping", &ValidatedIntegrationRequest::has_initial_mapping)
        .def("initial_mapping", &ValidatedIntegrationRequest::initial_mapping, py::return_value_policy::copy);

    py::class_<InputApproximationInterface>(module, "InputApproximationInterface")
        .def("is_reach_rigorous", &InputApproximationInterface::is_reach_rigorous);
    py::class_<ZeroApproximation, InputApproximationInterface>(module, "ZeroApproximation").def(py::init<>());
    py::class_<ConstantApproximation, InputApproximationInterface>(module, "ConstantApproximation").def(py::init<>());
    py::class_<AffineApproximation, InputApproximationInterface>(module, "AffineApproximation").def(py::init<>());
    py::class_<SinusoidalApproximation, InputApproximationInterface>(module, "SinusoidalApproximation").def(py::init<>());
    py::class_<PiecewiseApproximation, InputApproximationInterface>(module, "PiecewiseApproximation").def(py::init<>());
    py::class_<InputApproximation>(module, "InputApproximation")
        .def(py::init([](const InputApproximationInterface& approximation) { return InputApproximation(approximation.clone()); }))
        .def("is_reach_rigorous", &InputApproximation::is_reach_rigorous);

    py::class_<MultifunctionIntegratorInterface>(module, "MultifunctionIntegratorInterface");
    py::class_<TimedMultifunctionIntegratorInterface>(module, "TimedMultifunctionIntegratorInterface");
    py::class_<FirstOrderMultifunctionIntegrator, MultifunctionIntegratorInterface, TimedMultifunctionIntegratorInterface>(module, "FirstOrderMultifunctionIntegrator")
        .def(py::init<>());
    py::class_<SelectedParametricMultifunctionIntegrator, MultifunctionIntegratorInterface>(module, "SelectedParametricMultifunctionIntegrator")
        .def(py::init<const IntegratorBase&, const InputApproximation&, const InputApproximation&>(),
             py::arg("integrator"), py::arg("reach_approximation"), py::arg("evolve_approximation"))
        .def("reach_index", &SelectedParametricMultifunctionIntegrator::reach_index)
        .def("evolve_index", &SelectedParametricMultifunctionIntegrator::evolve_index);

    py::class_<ValidatedStepIntegratorInterface>(module, "ValidatedStepIntegratorInterface")
        .def("step", &ValidatedStepIntegratorInterface::step);
    py::class_<FlowStepIntegratorAdapter, ValidatedStepIntegratorInterface>(module, "FlowStepIntegratorAdapter")
        .def(py::init<const IntegratorInterface&>(), py::keep_alive<1, 2>());
    py::class_<MultiflowStepIntegratorAdapter, ValidatedStepIntegratorInterface>(module, "MultiflowStepIntegratorAdapter")
        .def(py::init<const MultifunctionIntegratorInterface&>(), py::keep_alive<1, 2>());
    py::class_<StructuredInclusionStepIntegratorAdapter, ValidatedStepIntegratorInterface>(module, "StructuredInclusionStepIntegratorAdapter")
        .def(py::init<const SelectedParametricMultifunctionIntegrator&>(), py::keep_alive<1, 2>());
}
