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
#include <stdexcept>
#include <utility>

#include "geometry/box.hpp"
#include "integration/integrator_interface.hpp"

namespace Ariadne {

/** @brief Common, validated enclosures of a single integration step.
 *
 *  reach() contains every intermediate state; evolve() contains every final
 *  state. These two independently certified bounds must not be conflated.
 *
 *  Algorithms which compute richer validated flow or inclusion patches retain
 *  that entire original result in a type-erased representation. Callers that
 *  need its parametric structure can retrieve it through representation_as<T>();
 *  no projection to boxes is forced on the underlying integration algorithm.
 *
 *  This is a numerical result, not a representation of a dynamical system.
 *  It deliberately does not choose or retry an integration method.
 */
class ValidatedIntegrationStep {
    StepSizeType _time_step;
    UpperBoxType _reach;
    UpperBoxType _evolve;
    std::any _representation;

  public:
    ValidatedIntegrationStep(const StepSizeType& time_step,
                             const UpperBoxType& reach,
                             const UpperBoxType& evolve,
                             std::any representation = {})
        : _time_step(time_step), _reach(reach), _evolve(evolve),
          _representation(std::move(representation)) {
        if (not (time_step > StepSizeType(0u)) ||
            reach.dimension() != evolve.dimension()) {
            throw std::invalid_argument("Validated integration step has incompatible bounds or duration");
        }
    }

    const StepSizeType& time_step() const { return _time_step; }
    const UpperBoxType& reach() const { return _reach; }
    const UpperBoxType& evolve() const { return _evolve; }

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
