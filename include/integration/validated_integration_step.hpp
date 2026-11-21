/***************************************************************************
 *            integration/validated_integration_step.hpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file integration/validated_integration_step.hpp
 *  @brief Common validated step result independent of the numerical method.
 */

#ifndef ARIADNE_VALIDATED_INTEGRATION_STEP_HPP
#define ARIADNE_VALIDATED_INTEGRATION_STEP_HPP

#include <any>
#include <optional>
#include <stdexcept>
#include <utility>

#include "function/function_patch.hpp"
#include "geometry/box.hpp"
#include "integration/integrator_interface.hpp"

namespace Ariadne {

/** @brief Common, validated enclosures of a single integration step.
 *
 *  reach() contains every intermediate state; evolve() contains every final
 *  state. These two independently certified bounds must not be conflated.
 *
 *  Algorithms which compute richer validated flow or inclusion patches expose
 *  their certified reach/endpoint maps through the same typed interface.
 *  Reach patches certify a tube; an endpoint patch certifies only the end.
 *  The original algorithm-specific result remains optional provenance and
 *  is not required for generic propagation or reconditioning.
 *
 *  This is a numerical result, not a representation of a dynamical system.
 *  It deliberately does not choose or retry an integration method.
 */
class ValidatedIntegrationStep {
    StepSizeType _time_step;
    StepSizeType _initial_time;
    UpperBoxType _reach;
    UpperBoxType _evolve;
    List<ValidatedVectorMultivariateFunctionPatch> _reach_models;
    std::optional<ValidatedVectorMultivariateFunctionPatch> _evolve_model;
    std::optional<SizeType> _reach_time_index;
    StepSizeType _reach_time_offset;
    bool _has_explicit_time;
    ExactBoxType _constant_parameter_domain;
    std::any _representation;

  public:
    ValidatedIntegrationStep(const StepSizeType& time_step,
                             const UpperBoxType& reach,
                             const UpperBoxType& evolve,
                             std::any representation = {},
                             List<ValidatedVectorMultivariateFunctionPatch> reach_models = {},
                             std::optional<ValidatedVectorMultivariateFunctionPatch> evolve_model = std::nullopt,
                             const StepSizeType& initial_time = StepSizeType(0u),
                             std::optional<SizeType> reach_time_index = std::nullopt,
                             bool has_explicit_time = false,
                             ExactBoxType constant_parameter_domain = ExactBoxType(0u),
                             const StepSizeType& reach_time_offset = StepSizeType(0u));

    const StepSizeType& time_step() const;
    const StepSizeType& initial_time() const;
    StepSizeType final_time() const;
    /** @brief Whether the numerical field uses an explicit time coordinate.
     *  This is numerical continuation metadata, not a system classification.
     */
    bool has_explicit_time() const;
    /** @brief Trajectory-constant parameters shared across integration steps.
     *  Their coordinates occupy the prefix immediately after original states.
     *  Later inclusion steps may append independent time-varying inputs.
     *  Structured maps must preserve this domain after removing reach time;
     *  inconsistent metadata is rejected when the step is constructed.
     */
    const ExactBoxType& constant_parameter_domain() const;
    const UpperBoxType& reach() const;
    const UpperBoxType& evolve() const;

    /** @brief Structured reach-tube certificates, if the algorithm provides any.
     *  The list may contain several patches with different parameter domains.
     *  An empty list means the algorithm certified only box enclosures.
     */
    const List<ValidatedVectorMultivariateFunctionPatch>& reach_models() const;
    bool has_reach_models() const;
    /** @brief Time-coordinate index in each reach patch, if maps are provided.
     *  The time coordinate may differ from the state dimension when a prior
     *  structured mapping has additional uncertainty parameters.
     */
    std::optional<SizeType> reach_time_index() const;

    /** @brief Add this numerical offset to patch time to obtain absolute time.
     *  Autonomous flow solvers can use step-local time; nonautonomous
     *  flows and inclusion solvers ordinarily use absolute time.
     */
    const StepSizeType& reach_time_offset() const;

    /** @brief Validated state sections at an absolute instant inside this step.
     *
     *  Returns one state patch for each reach patch containing that instant.
     *  Retains the original state and input parameter coordinates and removes
     *  only the numerical time argument. No box restart or solver-specific
     *  downcast is performed. An empty reach structure is not fabricated.
     */
    List<ValidatedVectorMultivariateFunctionPatch>
        reach_sections(const StepSizeType& absolute_time) const;

    /** @brief Endpoint mapping suitable for structured continuation.
     *  An integrator computing only box enclosures does not fabricate a map.
     *  Callers may use has_evolve_model() to choose a sound box-based path.
     */
    bool has_evolve_model() const;
    const ValidatedVectorMultivariateFunctionPatch& evolve_model() const;

    bool has_representation() const;

    /** @brief Recover an algorithm-specific validated result, when available.
     *  The returned pointer is owned by this step and remains valid for its
     *  lifetime. A mismatched type produces nullptr, never a lossy conversion.
     */
    template<class T>
    const T* representation_as() const {
        return std::any_cast<T>(&_representation);
    }
};

} // namespace Ariadne

#endif /* ARIADNE_VALIDATED_INTEGRATION_STEP_HPP */
