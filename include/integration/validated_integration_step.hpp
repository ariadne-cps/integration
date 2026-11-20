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
    UpperBoxType _reach;
    UpperBoxType _evolve;
    List<ValidatedVectorMultivariateFunctionPatch> _reach_models;
    std::optional<ValidatedVectorMultivariateFunctionPatch> _evolve_model;
    std::any _representation;

  public:
    ValidatedIntegrationStep(const StepSizeType& time_step,
                             const UpperBoxType& reach,
                             const UpperBoxType& evolve,
                             std::any representation = {},
                             List<ValidatedVectorMultivariateFunctionPatch> reach_models = {},
                             std::optional<ValidatedVectorMultivariateFunctionPatch> evolve_model = std::nullopt)
        : _time_step(time_step), _reach(reach), _evolve(evolve),
          _reach_models(std::move(reach_models)),
          _evolve_model(std::move(evolve_model)),
          _representation(std::move(representation)) {
        if (not (time_step > StepSizeType(0u)) ||
            reach.dimension() != evolve.dimension()) {
            throw std::invalid_argument("Validated integration step has incompatible bounds or duration");
        }
        for (const auto& model: _reach_models) {
            if (model.result_size() != reach.dimension()) {
                throw std::invalid_argument("Validated reach model has an incompatible state dimension");
            }
        }
        if (_evolve_model && _evolve_model->result_size() != evolve.dimension()) {
            throw std::invalid_argument("Validated endpoint model has an incompatible state dimension");
        }
    }

    const StepSizeType& time_step() const { return _time_step; }
    const UpperBoxType& reach() const { return _reach; }
    const UpperBoxType& evolve() const { return _evolve; }

    /** @brief Structured reach-tube certificates, if the algorithm provides any.
     *  The list may contain several patches with different parameter domains.
     *  An empty list means the algorithm certified only box enclosures.
     */
    const List<ValidatedVectorMultivariateFunctionPatch>& reach_models() const {
        return _reach_models;
    }
    bool has_reach_models() const { return not _reach_models.empty(); }

    /** @brief Endpoint mapping suitable for structured continuation.
     *  An integrator computing only box enclosures does not fabricate a map.
     *  Callers may use has_evolve_model() to choose a sound box-based path.
     */
    bool has_evolve_model() const { return _evolve_model.has_value(); }
    const ValidatedVectorMultivariateFunctionPatch& evolve_model() const {
        if (not _evolve_model) {
            throw std::logic_error("This integration step has no structured endpoint model");
        }
        return *_evolve_model;
    }

    bool has_representation() const { return _representation.has_value(); }

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
