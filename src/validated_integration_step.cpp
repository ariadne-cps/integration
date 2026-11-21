/***************************************************************************
 *            src/validated_integration_step.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/validated_integration_step.cpp
 *  @brief Implementation of validated integration step results.
 */

#include "integration/validated_integration_step.hpp"

namespace Ariadne {

ValidatedIntegrationStep::ValidatedIntegrationStep(const StepSizeType& time_step,
                                                   const UpperBoxType& reach,
                                                   const UpperBoxType& evolve,
                                                   std::any representation,
                                                   List<ValidatedVectorMultivariateFunctionPatch> reach_models,
                                                   std::optional<ValidatedVectorMultivariateFunctionPatch> evolve_model,
                                                   const StepSizeType& initial_time,
                                                   std::optional<SizeType> reach_time_index,
                                                   bool has_explicit_time,
                                                   ExactBoxType constant_parameter_domain,
                                                   const StepSizeType& reach_time_offset) : _time_step(time_step), _initial_time(initial_time), _reach(reach), _evolve(evolve),
                                                                                            _reach_models(std::move(reach_models)),
                                                                                            _evolve_model(std::move(evolve_model)),
                                                                                            _reach_time_index(reach_time_index ? reach_time_index
                                                                                                                               : (_reach_models.empty() ? std::nullopt : std::optional<SizeType>(reach.dimension()))),
                                                                                            _reach_time_offset(reach_time_offset),
                                                                                            _has_explicit_time(has_explicit_time),
                                                                                            _constant_parameter_domain(std::move(constant_parameter_domain)),
                                                                                            _representation(std::move(representation)) {
    if (not(time_step > StepSizeType(0u)) ||
        reach.dimension() != evolve.dimension()) {
        throw std::invalid_argument("Validated integration step has incompatible bounds or duration");
    }
    const SizeType n = reach.dimension();
    const SizeType np = _constant_parameter_domain.dimension();
    // Fixed parameters occupy the prefix after the original state
    // coordinates in endpoint maps and in reach maps AFTER removing time.
    // Input-approximation coordinates may differ between reach and evolve.
    auto validate_parameters = [&](const ValidatedVectorMultivariateFunctionPatch& model,
                                   std::optional<SizeType> time_index) {
        if (np == 0u) {
            return;
        }
        const SizeType nt = time_index ? 1u : 0u;
        if (model.argument_size() < n + np + nt || (time_index && *time_index < n)) {
            throw std::invalid_argument("Structured integration map has too few state or fixed-parameter coordinates");
        }
        const auto domain = model.domain();
        for (SizeType j = 0u; j != np; ++j) {
            const SizeType coordinate = n + j;
            const SizeType index = coordinate + (time_index && *time_index <= coordinate ? 1u : 0u);
            if (not(domain[index] == _constant_parameter_domain[j])) {
                throw std::invalid_argument("Structured integration map does not preserve the fixed-parameter domain");
            }
        }
    };
    for (const auto& model: _reach_models) {
        if (model.result_size() != reach.dimension() ||
            not _reach_time_index || *_reach_time_index >= model.argument_size()) {
            throw std::invalid_argument("Validated reach model has incompatible state or time coordinates");
        }
        validate_parameters(model, _reach_time_index);
        // Backend time may be local to the step or already absolute.
        // Each certified patch must lie within the correctly shifted
        // step interval; different patches may cover different parts.
        const IntervalDomainType model_step_time(
            _initial_time - _reach_time_offset,
            _initial_time + _time_step - _reach_time_offset);
        if (not definitely(subset(model.domain()[*_reach_time_index], model_step_time))) {
            throw std::invalid_argument("Validated reach patch lies outside the step's time interval");
        }
        // The common reach box must enclose ALL states certified by
        // every reach patch. A contradictory box is not a valid
        // lossy projection of a structured reach certificate.
        if (not definitely(subset(cast_exact_box(UpperBoxType(model.range())),
                                  cast_exact_box(_reach)))) {
            throw std::invalid_argument("Validated reach box does not enclose a reach patch");
        }
    }
    if (not _reach_models.empty()) {
        // Individual patch intervals may be split, overlap and arrive
        // unsorted. The structured reach must nonetheless certify a
        // section throughout the WHOLE step. Otherwise callers could
        // silently miss trajectories at intermediate physical times.
        const IntervalDomainType model_step_time(
            _initial_time - _reach_time_offset,
            _initial_time + _time_step - _reach_time_offset);
        auto covered = model_step_time.lower_bound();
        while (covered < model_step_time.upper_bound()) {
            auto next = covered;
            for (const auto& model: _reach_models) {
                // domain() returns a box BY VALUE. Retaining a reference
                // to one of its elements after this expression would
                // dangle, corrupting the time coverage comparison.
                const auto domain = model.domain();
                const auto& interval = domain[*_reach_time_index];
                if (interval.lower_bound() <= covered && interval.upper_bound() > next) {
                    next = interval.upper_bound();
                }
            }
            if (not(next > covered)) {
                throw std::invalid_argument(
                    "Validated reach patches do not cover the whole step time interval");
            }
            covered = next;
        }
    }
    if (_evolve_model && _evolve_model->result_size() != evolve.dimension()) {
        throw std::invalid_argument("Validated endpoint model has an incompatible state dimension");
    }
    if (_evolve_model) {
        validate_parameters(*_evolve_model, std::nullopt);
    }
    // Reach and endpoint may use different approximations and
    // different parameter domains. Validate their boxes separately.
    if (_evolve_model &&
        not definitely(subset(cast_exact_box(UpperBoxType(_evolve_model->range())),
                              cast_exact_box(_evolve)))) {
        throw std::invalid_argument("Validated endpoint box does not enclose the endpoint model");
    }
}

const StepSizeType& ValidatedIntegrationStep::time_step() const { return _time_step; }

const StepSizeType& ValidatedIntegrationStep::initial_time() const { return _initial_time; }

StepSizeType ValidatedIntegrationStep::final_time() const { return _initial_time + _time_step; }

bool ValidatedIntegrationStep::has_explicit_time() const { return _has_explicit_time; }

const ExactBoxType& ValidatedIntegrationStep::constant_parameter_domain() const {
    return _constant_parameter_domain;
}

const UpperBoxType& ValidatedIntegrationStep::reach() const { return _reach; }

const UpperBoxType& ValidatedIntegrationStep::evolve() const { return _evolve; }

const List<ValidatedVectorMultivariateFunctionPatch>& ValidatedIntegrationStep::reach_models() const {
    return _reach_models;
}

bool ValidatedIntegrationStep::has_reach_models() const { return not _reach_models.empty(); }

std::optional<SizeType> ValidatedIntegrationStep::reach_time_index() const { return _reach_time_index; }

const StepSizeType& ValidatedIntegrationStep::reach_time_offset() const { return _reach_time_offset; }

List<ValidatedVectorMultivariateFunctionPatch> ValidatedIntegrationStep::reach_sections(const StepSizeType& absolute_time) const {
    if (_reach_models.empty()) {
        throw std::logic_error("Box-only integration has no structured reach sections");
    }
    if (absolute_time < _initial_time || absolute_time > this->final_time()) {
        throw std::out_of_range("Reach section time lies outside the integration step");
    }
    const StepSizeType model_time = absolute_time - _reach_time_offset;
    const IntervalDomainType point(model_time, model_time);
    List<ValidatedVectorMultivariateFunctionPatch> sections;
    for (const auto& model: _reach_models) {
        if (definitely(subset(point, model.domain()[*_reach_time_index]))) {
            sections.append(partial_evaluate(model, *_reach_time_index, model_time));
        }
    }
    if (sections.empty()) {
        throw std::out_of_range("No certified reach patch contains the requested time");
    }
    return sections;
}

bool ValidatedIntegrationStep::has_evolve_model() const { return _evolve_model.has_value(); }

const ValidatedVectorMultivariateFunctionPatch& ValidatedIntegrationStep::evolve_model() const {
    if (not _evolve_model) {
        throw std::logic_error("This integration step has no structured endpoint model");
    }
    return *_evolve_model;
}

bool ValidatedIntegrationStep::has_representation() const { return _representation.has_value(); }

} // namespace Ariadne
