/***************************************************************************
 *            integration/integration_step_adapters.hpp
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

/** @file integration/integration_step_adapters.hpp
 *  @brief Lossless bridges from existing numerical step results.
 */

#ifndef ARIADNE_INTEGRATION_STEP_ADAPTERS_HPP
#define ARIADNE_INTEGRATION_STEP_ADAPTERS_HPP

#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "integration/validated_integration_step.hpp"

namespace Ariadne {

/** @brief Preserve the ordinary flow patch while exposing reach and endpoint.
 *
 *  The time domain must match the certified patch's time coordinate, so
 *  callers cannot relabel an already-computed flow with an arbitrary duration.
 *  This works for nonautonomous flows with a nonzero absolute initial time.
 */
ValidatedIntegrationStep validated_integration_step(
    const FlowStepModelType& flow,
    const Interval<StepSizeType>& time_domain,
    const StepSizeType& initial_time = StepSizeType(0u),
    std::optional<SizeType> time_coordinate = std::nullopt,
    bool has_explicit_time = false,
    const ExactBoxType& constant_parameter_domain = ExactBoxType(0u));

/** @brief Preserve the full structured inclusion result, including all
 *  reach patches and the parametrized endpoint map.
 */
ValidatedIntegrationStep validated_integration_step(
    const SelectedInclusionStep& step,
    const StepSizeType& initial_time = StepSizeType(0u),
    bool has_explicit_time = false,
    const ExactBoxType& constant_parameter_domain = ExactBoxType(0u),
    std::optional<SizeType> time_coordinate = std::nullopt);

/** @brief Box-only algorithms still satisfy the common result contract. */
ValidatedIntegrationStep validated_integration_step(
    const ValidatedMultiflowStep& step,
    const StepSizeType& initial_time = StepSizeType(0u),
    bool has_explicit_time = false,
    const ExactBoxType& constant_parameter_domain = ExactBoxType(0u));

} // namespace Ariadne

#endif /* ARIADNE_INTEGRATION_STEP_ADAPTERS_HPP */
