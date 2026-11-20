/***************************************************************************
 *            src/integration_step_adapters.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/integration_step_adapters.cpp
 *  @brief Validated step adapters from numerical results.
 */

#include "integration/integration_step_adapters.hpp"

namespace Ariadne {

ValidatedIntegrationStep validated_integration_step(
    const FlowStepModelType& flow,
    const Interval<StepSizeType>& time_domain,
    const StepSizeType& initial_time,
    std::optional<SizeType> time_coordinate,
    bool has_explicit_time,
    const ExactBoxType& constant_parameter_domain) {
    const SizeType n = flow.result_size();
    const SizeType time_index = time_coordinate.value_or(n);
    if (has_explicit_time && initial_time != time_domain.lower_bound()) {
        throw std::invalid_argument("Explicit-time flow cannot be relabelled with a different absolute start");
    }
    if (flow.argument_size() <= time_index ||
        not (flow.domain()[time_index] ==
             IntervalDomainType(time_domain.lower_bound(), time_domain.upper_bound()))) {
        throw std::invalid_argument("Ordinary flow patch has an incompatible time domain");
    }
    const StepSizeType duration = time_domain.upper_bound() - time_domain.lower_bound();
    const auto endpoint = partial_evaluate(flow, time_index, time_domain.upper_bound());
    List<ValidatedVectorMultivariateFunctionPatch> reach_models;
    reach_models.append(flow);
    return ValidatedIntegrationStep(duration, UpperBoxType(flow.range()),
                                    UpperBoxType(endpoint.range()), std::any(flow),
                                    std::move(reach_models), endpoint, initial_time,
                                    time_index, has_explicit_time,
                                    constant_parameter_domain,
                                    initial_time - time_domain.lower_bound());
}

ValidatedIntegrationStep validated_integration_step(
    const SelectedInclusionStep& step,
    const StepSizeType& initial_time,
    bool has_explicit_time,
    const ExactBoxType& constant_parameter_domain,
    std::optional<SizeType> time_coordinate) {
    const auto bounds = step.bounding_step();
    return ValidatedIntegrationStep(bounds.time_step, bounds.reach,
                                    bounds.evolve, std::any(step),
                                    step.reach_models, step.evolve_model, initial_time,
                                    time_coordinate, has_explicit_time,
                                    constant_parameter_domain);
}

ValidatedIntegrationStep validated_integration_step(
    const ValidatedMultiflowStep& step,
    const StepSizeType& initial_time,
    bool has_explicit_time,
    const ExactBoxType& constant_parameter_domain) {
    return ValidatedIntegrationStep(step.time_step, step.reach, step.evolve,
                                    {}, {}, std::nullopt, initial_time,
                                    std::nullopt, has_explicit_time,
                                    constant_parameter_domain);
}

} // namespace Ariadne
