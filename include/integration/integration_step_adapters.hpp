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
inline ValidatedIntegrationStep validated_integration_step(
    const FlowStepModelType& flow,
    const Interval<StepSizeType>& time_domain,
    const StepSizeType& initial_time = StepSizeType(0u)) {
    const SizeType n = flow.result_size();
    if (flow.argument_size() < n + 1u ||
        not (flow.domain()[n] ==
             IntervalDomainType(time_domain.lower_bound(), time_domain.upper_bound()))) {
        throw std::invalid_argument("Ordinary flow patch has an incompatible time domain");
    }
    const StepSizeType duration = time_domain.upper_bound() - time_domain.lower_bound();
    const auto endpoint = partial_evaluate(flow, n, time_domain.upper_bound());
    List<ValidatedVectorMultivariateFunctionPatch> reach_models;
    reach_models.append(flow);
    return ValidatedIntegrationStep(duration, UpperBoxType(flow.range()),
                                    UpperBoxType(endpoint.range()), std::any(flow),
                                    std::move(reach_models), endpoint, initial_time);
}

/** @brief Preserve the full structured inclusion result, including all
 *  reach patches and the parametrized endpoint map.
 */
inline ValidatedIntegrationStep validated_integration_step(
    const SelectedInclusionStep& step,
    const StepSizeType& initial_time = StepSizeType(0u)) {
    const auto bounds = step.bounding_step();
    return ValidatedIntegrationStep(bounds.time_step, bounds.reach,
                                    bounds.evolve, std::any(step),
                                    step.reach_models, step.evolve_model, initial_time);
}

/** @brief Box-only algorithms still satisfy the common result contract. */
inline ValidatedIntegrationStep validated_integration_step(
    const ValidatedMultiflowStep& step,
    const StepSizeType& initial_time = StepSizeType(0u)) {
    return ValidatedIntegrationStep(step.time_step, step.reach, step.evolve,
                                    {}, {}, std::nullopt, initial_time);
}

} // namespace Ariadne

#endif /* ARIADNE_INTEGRATION_STEP_ADAPTERS_HPP */
