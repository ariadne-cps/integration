/***************************************************************************
 *            src/validated_step_integrator.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/validated_step_integrator.cpp
 *  @brief Validated numerical integration requests and solver bridges.
 */

#include "integration/validated_step_integrator.hpp"

namespace Ariadne {

ValidatedIntegrationRequest::ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                                         const ExactBoxType& initial,
                                                         const Suggestion<StepSizeType>& maximum_step,
                                                         const StepSizeType& initial_time)
    : _field(field), _initial(initial), _time_step(static_cast<const StepSizeType&>(maximum_step)),
      _initial_time(initial_time), _maximum_step(true) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                                         const ValidatedIntegrationStep& previous,
                                                         const Suggestion<StepSizeType>& maximum_step)
    : ValidatedIntegrationRequest(field, previous, static_cast<const StepSizeType&>(maximum_step), previous.evolve()) {
    _maximum_step = true;
    _bounding_box.reset();
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                                         const ExactBoxType& initial,
                                                         const StepSizeType& time_step,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(initial), _time_step(time_step),
                                                                                             _initial_time(0u), _bounding_box(bounding_box) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedBoxMultifunctionInterface> field,
                                                         const ExactBoxType& initial,
                                                         const StepSizeType& time_step,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(initial), _time_step(time_step),
                                                                                             _initial_time(0u), _bounding_box(bounding_box) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                                         const ExactBoxType& initial,
                                                         const StepSizeType& initial_time,
                                                         const StepSizeType& time_step,
                                                         const ExactBoxType& constant_parameters,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(initial), _time_step(time_step),
                                                                                             _initial_time(initial_time), _bounding_box(bounding_box),
                                                                                             _has_explicit_time(true),
                                                                                             _constant_parameter_domain(constant_parameters) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedTimedBoxMultifunctionInterface> field,
                                                         const ExactBoxType& initial,
                                                         const StepSizeType& initial_time,
                                                         const StepSizeType& time_step,
                                                         const ExactBoxType& constant_parameters,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(initial), _time_step(time_step),
                                                                                             _initial_time(initial_time), _bounding_box(bounding_box),
                                                                                             _has_explicit_time(true),
                                                                                             _constant_parameter_domain(constant_parameters) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                                         const ValidatedIntegrationStep& previous,
                                                         const StepSizeType& time_step,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(continuation_domain(previous)),
                                                                                             _time_step(time_step), _initial_time(previous.final_time()),
                                                                                             _bounding_box(bounding_box),
                                                                                             _initial_mapping(previous.has_evolve_model()
                                                                                                                  ? std::optional<ValidatedVectorMultivariateFunctionPatch>(previous.evolve_model())
                                                                                                                  : std::nullopt),
                                                                                             _has_explicit_time(previous.has_explicit_time()),
                                                                                             _constant_parameter_domain(previous.constant_parameter_domain()) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedBoxMultifunctionInterface> field,
                                                         const ValidatedIntegrationStep& previous,
                                                         const StepSizeType& time_step,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(continuation_domain(previous)),
                                                                                             _time_step(time_step), _initial_time(previous.final_time()),
                                                                                             _bounding_box(bounding_box),
                                                                                             _initial_mapping(previous.has_evolve_model()
                                                                                                                  ? std::optional<ValidatedVectorMultivariateFunctionPatch>(previous.evolve_model())
                                                                                                                  : std::nullopt),
                                                                                             _has_explicit_time(previous.has_explicit_time()),
                                                                                             _constant_parameter_domain(previous.constant_parameter_domain()) {
    this->validate();
}

ValidatedIntegrationRequest::ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedTimedBoxMultifunctionInterface> field,
                                                         const ValidatedIntegrationStep& previous,
                                                         const StepSizeType& time_step,
                                                         const UpperBoxType& bounding_box) : _field(field), _initial(continuation_domain(previous)),
                                                                                             _time_step(time_step), _initial_time(previous.final_time()),
                                                                                             _bounding_box(bounding_box),
                                                                                             _initial_mapping(previous.has_evolve_model()
                                                                                                                  ? std::optional<ValidatedVectorMultivariateFunctionPatch>(previous.evolve_model())
                                                                                                                  : std::nullopt),
                                                                                             _has_explicit_time(true),
                                                                                             _constant_parameter_domain(previous.constant_parameter_domain()) {
    this->validate();
}

ExactBoxType ValidatedIntegrationRequest::continuation_domain(const ValidatedIntegrationStep& previous) {
    if (previous.has_evolve_model()) {
        const SizeType n = previous.evolve().dimension();
        if (previous.evolve_model().argument_size() < n) {
            throw std::invalid_argument("Previous endpoint mapping has too few state coordinates");
        }
        return ExactBoxType(previous.evolve_model().domain()[range(0u, n)]);
    }
    return cast_exact_box(previous.evolve());
}

const ValidatedIntegrationRequest::Field& ValidatedIntegrationRequest::field() const { return _field; }

const ExactBoxType& ValidatedIntegrationRequest::initial() const { return _initial; }

const StepSizeType& ValidatedIntegrationRequest::time_step() const { return _time_step; }

const StepSizeType& ValidatedIntegrationRequest::initial_time() const { return _initial_time; }

bool ValidatedIntegrationRequest::step_size_is_maximum() const { return _maximum_step; }

bool ValidatedIntegrationRequest::has_bounding_box() const { return _bounding_box.has_value(); }

const UpperBoxType& ValidatedIntegrationRequest::bounding_box() const {
    if (!has_bounding_box()) {
        throw std::logic_error("The adaptive request leaves the state bound to the solver");
    }
    return *_bounding_box;
}

bool ValidatedIntegrationRequest::has_explicit_time() const { return _has_explicit_time; }

const ExactBoxType& ValidatedIntegrationRequest::constant_parameter_domain() const {
    return _constant_parameter_domain;
}

bool ValidatedIntegrationRequest::has_initial_mapping() const { return _initial_mapping.has_value(); }

const ValidatedVectorMultivariateFunctionPatch& ValidatedIntegrationRequest::initial_mapping() const {
    if (not _initial_mapping) {
        throw std::logic_error("This integration request has no initial mapping");
    }
    return *_initial_mapping;
}

void ValidatedIntegrationRequest::validate() const {
    if (_maximum_step && _has_explicit_time) {
        throw std::invalid_argument("Adaptive step requests require an autonomous equation field");
    }
    if (not(_time_step > StepSizeType(0u)) ||
        (_bounding_box && _initial.dimension() != _bounding_box->dimension())) {
        throw std::invalid_argument("Validated integration request has incompatible dimensions or duration");
    }
    if (_constant_parameter_domain.dimension() != 0u && not _has_explicit_time) {
        throw std::invalid_argument("Trajectory-constant equation parameters require explicit-time field semantics");
    }
    if (_initial_mapping &&
        (_initial_mapping->result_size() != _initial.dimension() ||
         _initial_mapping->argument_size() < _initial.dimension() ||
         not(_initial_mapping->domain()[range(0u, _initial.dimension())] == _initial))) {
        throw std::invalid_argument("Validated continuation mapping has an incompatible initial-state domain");
    }
    if (_initial_mapping && _constant_parameter_domain.dimension() != 0u) {
        const SizeType n = _initial.dimension();
        const SizeType p = _constant_parameter_domain.dimension();
        if (_initial_mapping->argument_size() < n + p ||
            not(_initial_mapping->domain()[range(n, n + p)] ==
                _constant_parameter_domain)) {
            throw std::invalid_argument("Continuation must reuse the original trajectory-constant parameter domain");
        }
    }
}

FlowStepIntegratorAdapter::FlowStepIntegratorAdapter(const IntegratorInterface& integrator) : _integrator(integrator) {}

MultiflowStepIntegratorAdapter::MultiflowStepIntegratorAdapter(const MultifunctionIntegratorInterface& integrator) : _integrator(integrator) {}

StructuredInclusionStepIntegratorAdapter::StructuredInclusionStepIntegratorAdapter(
    const SelectedParametricMultifunctionIntegrator& integrator) : _integrator(integrator) {}

ValidatedIntegrationStep FlowStepIntegratorAdapter::step(const ValidatedIntegrationRequest& request) const {
    const auto* field = request.field_as<ValidatedVectorMultivariateFunction>();
    if (field == nullptr) {
        throw std::invalid_argument("This flow integrator requires a validated single-valued vector field");
    }
    const SizeType n = request.initial().dimension();
    const SizeType np = request.constant_parameter_domain().dimension();
    const bool explicit_time = request.has_explicit_time();
    if (field->result_size() != n ||
        field->argument_size() != n + (explicit_time ? 1u + np : 0u)) {
        throw std::invalid_argument("Equation field does not match state, time and parameter arguments");
    }
    const ExactBoxType current = request.has_initial_mapping()
                                   ? cast_exact_box(UpperBoxType(request.initial_mapping().range()))
                                   : request.initial();
    if (request.has_initial_mapping() && request.has_bounding_box() &&
        not definitely(subset(current, cast_exact_box(request.bounding_box())))) {
        throw std::invalid_argument("Equation continuation endpoint exceeds the state bound");
    }
    const auto local = [&]() -> FlowStepModelType {
        if (request.step_size_is_maximum()) {
            auto flow = _integrator.flow_step(*field, current, suggest(request.time_step()));
            auto domain = flow.domain();
            const StepSizeType h = min(request.time_step(), StepSizeType(domain[n].upper_bound()));
            if (!(h > StepSizeType(0u)) || domain[n].lower_bound() > 0) {
                throw std::invalid_argument("Adaptive solver returned no positive flow interval starting at zero");
            }
            domain[n] = IntervalDomainType(StepSizeType(0u), h);
            return restriction(flow, domain);
        }
        if (explicit_time) {
            return _integrator.flow_step(*field, current,
                                         Interval<StepSizeType>(request.initial_time(), request.initial_time() + request.time_step()),
                                         request.constant_parameter_domain(), request.bounding_box());
        }
        return _integrator.flow_step(*field, current, request.time_step(), request.bounding_box());
    }();
    const Interval<StepSizeType> time_domain = request.step_size_is_maximum()
                                                 ? Interval<StepSizeType>(StepSizeType(0u), StepSizeType(local.domain()[n].upper_bound()))
                                                 : Interval<StepSizeType>(explicit_time ? request.initial_time() : StepSizeType(0u),
                                                                          (explicit_time ? request.initial_time() : StepSizeType(0u)) + request.time_step());
    if (!request.has_initial_mapping()) {
        return validated_integration_step(local, time_domain, request.initial_time(),
                                          std::nullopt, explicit_time, request.constant_parameter_domain());
    }

    // Preserve the previous map and every existing uncertainty coordinate.
    const auto& previous = request.initial_mapping();
    const auto embedded = embed(previous, IntervalDomainType(time_domain));
    const SizeType time_index = previous.argument_size();
    const auto time_function = factory(embedded).create(
        ValidatedScalarMultivariateFunction::coordinate(time_index + 1u, time_index));
    auto arguments = join(embedded, time_function);
    for (SizeType i = 0u; i != np; ++i) {
        // Fixed parameters are existing original-domain coordinates;
        // do NOT introduce an independent value at every step.
        const auto parameter = factory(embedded).create(
            ValidatedScalarMultivariateFunction::coordinate(
                time_index + 1u, n + i));
        arguments = join(arguments, parameter);
    }
    const FlowStepModelType composed(compose(local, arguments));
    return validated_integration_step(composed, time_domain,
                                      request.initial_time(), time_index,
                                      explicit_time, request.constant_parameter_domain());
}

ValidatedIntegrationStep MultiflowStepIntegratorAdapter::step(const ValidatedIntegrationRequest& request) const {
    if (request.step_size_is_maximum()) {
        throw std::invalid_argument("This inclusion adapter requires an exact step and a state bound");
    }
    if (request.has_initial_mapping()) {
        throw std::invalid_argument("Box-only integration cannot discard a structured endpoint map");
    }
    if (const auto* timed = request.field_as<ValidatedTimedBoxMultifunctionInterface>()) {
        const auto* solver = dynamic_cast<const TimedMultifunctionIntegratorInterface*>(&_integrator);
        if (solver == nullptr) {
            throw std::invalid_argument("The configured inclusion integrator lacks timed-field capability");
        }
        return validated_integration_step(solver->flow_step(
                                              *timed, request.initial(), request.initial_time(), request.time_step(),
                                              request.constant_parameter_domain(), request.bounding_box()),
                                          request.initial_time(), true, request.constant_parameter_domain());
    }
    const auto* field = request.field_as<ValidatedBoxMultifunctionInterface>();
    if (field == nullptr) {
        throw std::invalid_argument("This inclusion integrator requires a validated box multifunction");
    }
    if (request.constant_parameter_domain().dimension() != 0u) {
        throw std::invalid_argument("Autonomous box field cannot consume fixed equation parameters");
    }
    return validated_integration_step(_integrator.flow_step(
                                          *field, request.initial(), request.time_step(), request.bounding_box()),
                                      request.initial_time(), request.has_explicit_time(),
                                      request.constant_parameter_domain());
}

ValidatedIntegrationStep StructuredInclusionStepIntegratorAdapter::step(const ValidatedIntegrationRequest& request) const {
    if (request.step_size_is_maximum()) {
        throw std::invalid_argument("This inclusion adapter requires an exact step and a state bound");
    }
    if (const auto* timed = request.field_as<ValidatedTimedBoxMultifunctionInterface>()) {
        const auto* previous = request.has_initial_mapping()
                                 ? &request.initial_mapping()
                                 : nullptr;
        const auto native = _integrator.timed_structured_step(
            *timed, request.initial(), request.initial_time(), request.time_step(),
            request.constant_parameter_domain(), request.bounding_box(), previous);
        const SizeType time_index =
            request.initial().dimension() + request.constant_parameter_domain().dimension();
        return validated_integration_step(
            native, request.initial_time(), true,
            request.constant_parameter_domain(), time_index);
    }
    const auto* field = request.field_as<ValidatedBoxMultifunctionInterface>();
    if (field == nullptr) {
        throw std::invalid_argument("Structured inclusion integration requires a validated box multifunction");
    }
    if (request.has_initial_mapping()) {
        return validated_integration_step(_integrator.structured_step(
                                              *field, request.initial(), request.initial_mapping(),
                                              request.initial_time(), request.time_step(), request.bounding_box()),
                                          request.initial_time(), request.has_explicit_time(),
                                          request.constant_parameter_domain());
    }
    return validated_integration_step(_integrator.structured_step(
                                          *field, request.initial(), request.initial_time(),
                                          request.time_step(), request.bounding_box()),
                                      request.initial_time(), request.has_explicit_time(),
                                      request.constant_parameter_domain());
}

} // namespace Ariadne
