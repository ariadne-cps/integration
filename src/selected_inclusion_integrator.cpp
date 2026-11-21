/***************************************************************************
 *            src/selected_inclusion_integrator.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/selected_inclusion_integrator.cpp
 *  @brief Selected high-order inclusion method implementations.
 */

#include "integration/selected_inclusion_integrator.hpp"

namespace Ariadne {

ValidatedMultiflowStep SelectedInclusionStep::bounding_step() const {
    if (reach_models.empty()) {
        throw FlowBoundsException("The selected reach method did not return a reach model");
    }
    UpperBoxType reach_box(reach_models[0u].range());
    for (SizeType i = 1u; i != reach_models.size(); ++i) {
        reach_box = hull(reach_box, UpperBoxType(reach_models[i].range()));
    }
    return ValidatedMultiflowStep{time_step, reach_box, UpperBoxType(evolve_model.range())};
}

SelectedInclusionIntegrator::SelectedInclusionIntegrator(const InclusionIntegratorFactory& factory,
                                                         const EffectiveVectorMultivariateFunction& field,
                                                         const BoxDomainType& inputs,
                                                         const InputApproximation& reach_approximation,
                                                         const InputApproximation& evolve_approximation)
    : _field(field),
      _inputs(inputs),
      _factory(factory),
      _reach_approximation(reach_approximation),
      _evolve_approximation(evolve_approximation) {
    if (field.argument_size() != field.result_size() + inputs.dimension()) {
        throw std::invalid_argument("Parametric inclusion field and input dimensions do not match");
    }
    if (not _reach_approximation.is_reach_rigorous()) {
        throw std::invalid_argument("The selected reach approximation does not certify the whole tube");
    }
}

Nat SelectedInclusionIntegrator::reach_index() const { return _reach_approximation.reference().index(); }

Nat SelectedInclusionIntegrator::evolve_index() const { return _evolve_approximation.reference().index(); }

SelectedInclusionStep SelectedInclusionIntegrator::step(const BoxDomainType& state_domain,
                                                        const ValidatedVectorMultivariateFunctionPatch& initial_mapping,
                                                        const UpperBoxType& state_bound,
                                                        const TimeStepType& initial_time,
                                                        const StepSizeType& time_step) const {
    if (not(time_step > StepSizeType(0u))) {
        throw std::invalid_argument("Inclusion integration requires a positive time step");
    }
    const SizeType n = _field.result_size();
    if (state_domain.dimension() != n || state_bound.dimension() != n ||
        initial_mapping.result_size() != n || initial_mapping.argument_size() < n) {
        throw std::invalid_argument("Inclusion step state and initial mapping dimensions do not match");
    }
    // Coordinates 0..n-1 encode the original state domain. Additional
    // arguments encode uncertainties retained from earlier steps.
    // The caller may append new parameters but must not rebase these.
    if (not(initial_mapping.domain()[range(0u, n)] == state_domain)) {
        throw std::invalid_argument("Inclusion step initial mapping has an incompatible state domain");
    }
    // The patch domain encodes original state coordinates and all
    // uncertainty parameters; it is NOT the current state enclosure.
    // The local flow must instead be certified on the image of the
    // incoming endpoint mapping. In a continuation these differ, and
    // reusing the original state box would evaluate the flow outside
    // its validated domain (or incorrectly discard the previous step).
    const ExactBoxType current_state = cast_exact_box(
        UpperBoxType(initial_mapping.range()));
    if (not definitely(subset(current_state, cast_exact_box(state_bound)))) {
        throw std::invalid_argument(
            "Current structured inclusion endpoint exceeds the certified state bound");
    }

    auto reach_method = _factory.create(_field, _inputs, _reach_approximation);
    if (not reach_method.is_reach_rigorous()) {
        throw std::invalid_argument("The selected reach implementation does not certify the whole tube");
    }
    auto reaches = reach_method.reach(current_state, initial_mapping, state_bound, initial_time, time_step);
    if (reaches.empty()) {
        throw FlowBoundsException("The selected reach implementation returned no reach model");
    }

    auto evolve_method = _factory.create(_field, _inputs, _evolve_approximation);
    const TimeStepType final_time = lower_bound(initial_time + time_step);
    if (reach_method.index() == evolve_method.index()) {
        if (reaches.size() != 1u) {
            throw FlowBoundsException("Endpoint evolution requires one model from the selected method");
        }
        auto endpoint = evolve_method.evolve(reaches[0u], final_time);
        return SelectedInclusionStep{time_step, std::move(reaches), std::move(endpoint)};
    }

    auto evolve_reaches = evolve_method.reach(current_state, initial_mapping, state_bound, initial_time, time_step);
    if (evolve_reaches.size() != 1u) {
        throw FlowBoundsException("The selected evolve method did not return exactly one endpoint source");
    }
    auto endpoint = evolve_method.evolve(evolve_reaches[0u], final_time);
    return SelectedInclusionStep{time_step, std::move(reaches), std::move(endpoint)};
}

ValidatedInputFieldBoxAdapterDP::ValidatedInputFieldBoxAdapterDP(const EffectiveVectorMultivariateFunction& field,
                                                                 const BoxDomainType& inputs,
                                                                 const BoxDomainType& state_domain)
    : _field(field), _inputs(inputs), _state_domain(state_domain) {
    if (inputs.dimension() == 0u || field.result_size() != state_domain.dimension() ||
        field.argument_size() != field.result_size() + inputs.dimension()) {
        throw std::invalid_argument("Parametric multifunction requires n states and at least one input");
    }
}

SizeType ValidatedInputFieldBoxAdapterDP::argument_size() const { return _field.result_size(); }

SizeType ValidatedInputFieldBoxAdapterDP::result_size() const { return _field.result_size(); }

BoxDomainType ValidatedInputFieldBoxAdapterDP::domain() const { return _state_domain; }

const EffectiveVectorMultivariateFunction& ValidatedInputFieldBoxAdapterDP::vector_field() const { return _field; }

const BoxDomainType& ValidatedInputFieldBoxAdapterDP::inputs() const { return _inputs; }

InclusionIntegratorFieldData ValidatedInputFieldBoxAdapterDP::inclusion_integrator_field() const {
    return InclusionIntegratorFieldData{_field, _inputs};
}

UpperBoxType ValidatedInputFieldBoxAdapterDP::evaluate(const UpperBoxType& box) const {
    if (box.dimension() != this->argument_size() ||
        not definitely(subset(box, _state_domain))) {
        throw std::invalid_argument("Parametric multifunction evaluation outside the state domain");
    }
    const SizeType n = this->argument_size();
    Vector<FloatDPBounds> arguments(n + _inputs.dimension(), [&](SizeType i) {
        if (i < n) {
            return FloatDPBounds(box[i].lower_bound(), box[i].upper_bound());
        }
        const auto& input = _inputs[i - n];
        return FloatDPBounds(input.lower_bound(), input.upper_bound());
    });
    const Vector<FloatDPBounds> values = _field(arguments);
    return UpperBoxType(values.size(), [&](SizeType i) {
        return UpperIntervalType(values[i].lower(), values[i].upper());
    });
}

SelectedParametricMultifunctionIntegrator::SelectedParametricMultifunctionIntegrator(const IntegratorBase& integrator,
                                                                                     const InputApproximation& reach_approximation,
                                                                                     const InputApproximation& evolve_approximation)
    : _factory(integrator),
      _function_factory(integrator.function_factory()),
      _reach_approximation(reach_approximation),
      _evolve_approximation(evolve_approximation) {
    if (not _reach_approximation.is_reach_rigorous()) {
        throw std::invalid_argument("High-order reach requires a tube-rigorous input approximation");
    }
}

Nat SelectedParametricMultifunctionIntegrator::reach_index() const { return _reach_approximation.reference().index(); }

Nat SelectedParametricMultifunctionIntegrator::evolve_index() const { return _evolve_approximation.reference().index(); }

SelectedInclusionStep SelectedParametricMultifunctionIntegrator::structured_step(const ValidatedBoxMultifunctionInterface& field,
                                                                                 const ExactBoxType& initial,
                                                                                 const ValidatedVectorMultivariateFunctionPatch& initial_mapping,
                                                                                 const TimeStepType& initial_time,
                                                                                 const StepSizeType& time_step,
                                                                                 const UpperBoxType& state_bound) const {
    if (initial.dimension() != field.argument_size() ||
        state_bound.dimension() != field.argument_size()) {
        throw std::invalid_argument("Parametric multifunction state dimensions do not match");
    }
    if (not definitely(subset(initial, field.domain())) ||
        not definitely(subset(state_bound, field.domain()))) {
        throw std::invalid_argument("High-order integration requires the initial and bounding boxes inside the state domain");
    }

    // Query a separate numerical capability, never a concrete field type.
    // A box-only field remains valid for other algorithms but cannot be
    // silently converted or dispatched to a different solver.
    const auto* provider = dynamic_cast<const InclusionIntegratorFieldProvider*>(&field);
    if (provider == nullptr) {
        throw std::invalid_argument("This multifunction does not provide high-order inclusion solver field data");
    }
    const auto form = provider->inclusion_integrator_field();
    if (form.field.result_size() != initial.dimension() ||
        form.field.argument_size() != initial.dimension() + form.inputs.dimension() ||
        form.inputs.dimension() == 0u) {
        throw std::invalid_argument("The inclusion solver field data has incompatible dimensions");
    }

    // The original initial domain is retained even when the current
    // endpoint map has additional independent input parameters.
    if (initial_mapping.result_size() != initial.dimension() ||
        initial_mapping.argument_size() < initial.dimension() ||
        not(initial_mapping.domain()[range(0u, initial.dimension())] == initial)) {
        throw std::invalid_argument("Structured continuation has incompatible mapping coordinates");
    }
    if (not definitely(subset(UpperBoxType(initial_mapping.range()),
                              cast_exact_box(state_bound)))) {
        throw std::invalid_argument("Structured continuation endpoint is outside the certified state bound");
    }
    // The legacy backend analyses input additivity on a formula
    // function. A vector assembled from scalar functions is equally
    // valid input, but its wrapper does not expose that concrete type.
    // Materialise its symbolic representation here, at the solver boundary.
    const Vector<EffectiveFormula> coordinates(form.field.argument_size(), [&](SizeType i) {
        return EffectiveFormula::coordinate(i);
    });
    const EffectiveVectorMultivariateFunction formula_field =
        EffectiveVectorMultivariateFormulaFunction(form.field.argument_size(), form.field(coordinates));
    SelectedInclusionIntegrator selected(_factory, formula_field, form.inputs,
                                         _reach_approximation, _evolve_approximation);
    return selected.step(initial, initial_mapping, state_bound, initial_time, time_step);
}

SelectedInclusionStep SelectedParametricMultifunctionIntegrator::structured_step(const ValidatedBoxMultifunctionInterface& field,
                                                                                 const ExactBoxType& initial,
                                                                                 const TimeStepType& initial_time,
                                                                                 const StepSizeType& time_step,
                                                                                 const UpperBoxType& state_bound) const {
    const auto identity = _function_factory.create_identity(initial);
    return this->structured_step(field, initial, identity,
                                 initial_time, time_step, state_bound);
}

SelectedInclusionStep SelectedParametricMultifunctionIntegrator::structured_step(const ValidatedBoxMultifunctionInterface& field,
                                                                                 const ExactBoxType& initial,
                                                                                 const StepSizeType& time_step,
                                                                                 const UpperBoxType& state_bound) const {
    return this->structured_step(field, initial, TimeStepType(0u),
                                 time_step, state_bound);
}

SelectedInclusionStep SelectedParametricMultifunctionIntegrator::timed_structured_step(
    const ValidatedTimedBoxMultifunctionInterface& field,
    const ExactBoxType& initial,
    const TimeStepType& initial_time,
    const StepSizeType& time_step,
    const ExactBoxType& parameters,
    const UpperBoxType& state_bound,
    const ValidatedVectorMultivariateFunctionPatch* previous) const {
    const SizeType n = initial.dimension();
    const SizeType np = parameters.dimension();
    const SizeType ns = n + 1u + np;
    const IntervalDomainType times(initial_time, initial_time + time_step);
    if (not(time_step > StepSizeType(0u)) || field.argument_size() != n ||
        field.result_size() != n || state_bound.dimension() != n ||
        field.parameter_domain().dimension() != np ||
        not definitely(subset(initial, field.domain())) ||
        not definitely(subset(state_bound, field.domain())) ||
        not definitely(subset(times, field.time_domain())) ||
        not definitely(subset(parameters, field.parameter_domain()))) {
        throw std::invalid_argument("Timed high-order inclusion request outside numerical field domain");
    }
    const auto* provider = dynamic_cast<const TimedInclusionIntegratorFieldProvider*>(&field);
    if (provider == nullptr) {
        throw std::invalid_argument("Timed field does not supply structured high-order solver data");
    }
    const auto data = provider->timed_inclusion_integrator_field();
    const SizeType args = ns + data.inputs.dimension();
    if (data.inputs.dimension() == 0u || data.field.result_size() != n ||
        data.field.argument_size() != args) {
        throw std::invalid_argument("Timed structured field has incompatible state, time or input dimensions");
    }

    // The legacy inclusion solver analyses input additivity directly
    // on EffectiveVectorMultivariateFormulaFunction. Building the
    // augmented field with join() erases that representation, even when
    // each original component is symbolic. Evaluate the supplied
    // structured field over symbolic coordinates and construct one
    // formula function for all state derivatives, including clock'=1
    // and parameter'=0. Only actual input coordinates are solver inputs.
    const Vector<EffectiveFormula> coordinates(args, [&](SizeType i) {
        return EffectiveFormula::coordinate(i);
    });
    const Vector<EffectiveFormula> original = data.field(coordinates);
    const Vector<EffectiveFormula> derivatives(ns, [&](SizeType i) {
        return i < n ? original[i] : (i == n ? EffectiveFormula::constant(EffectiveNumber(1u)) : EffectiveFormula::zero());
    });
    const EffectiveVectorMultivariateFunction lifted_field =
        EffectiveVectorMultivariateFormulaFunction(args, derivatives);

    const BoxDomainType clock_start(1u,
                                    IntervalDomainType(initial_time, initial_time));
    const BoxDomainType augmented_initial = product(
        BoxDomainType(initial), clock_start, BoxDomainType(parameters));
    const UpperBoxType clock_bound(1u, UpperIntervalType(
                                           field.time_domain().lower_bound(), field.time_domain().upper_bound()));
    const UpperBoxType augmented_bound = product(
        state_bound, clock_bound, UpperBoxType(field.parameter_domain()));

    ValidatedVectorMultivariateFunctionPatch initial_mapping;
    if (previous == nullptr) {
        initial_mapping = _function_factory.create_identity(augmented_initial);
    } else {
        const SizeType old_args = previous->argument_size();
        if (previous->result_size() != n || old_args < n + np ||
            not(previous->domain()[range(0u, n)] == BoxDomainType(initial)) ||
            not(previous->domain()[range(n, n + np)] == BoxDomainType(parameters))) {
            throw std::invalid_argument("Timed continuation must retain original state and parameter coordinates");
        }
        const BoxDomainType domain = product(
            BoxDomainType(initial), clock_start,
            previous->domain()[range(n, old_args)]);
        auto remap = _function_factory.create_zeros(old_args, domain);
        for (SizeType i = 0u; i != old_args; ++i) {
            remap[i] = _function_factory.create_coordinate(
                domain, i < n ? i : i + 1u);
        }
        initial_mapping = join(
            compose(*previous, remap),
            _function_factory.create_coordinate(domain, n));
        for (SizeType i = 0u; i != np; ++i) {
            initial_mapping = join(initial_mapping,
                                   _function_factory.create_coordinate(domain, n + 1u + i));
        }
    }

    SelectedInclusionIntegrator solver(_factory, lifted_field, data.inputs,
                                       _reach_approximation, _evolve_approximation);
    auto result = solver.step(
        augmented_initial, initial_mapping, augmented_bound, initial_time, time_step);
    List<ValidatedVectorMultivariateFunctionPatch> reach_models;
    for (const auto& model: result.reach_models) {
        auto copy = model;
        auto state = copy[range(0u, n)];
        reach_models.append(partial_evaluate(state, n, initial_time));
    }
    auto endpoint = result.evolve_model;
    auto state_endpoint = endpoint[range(0u, n)];
    return SelectedInclusionStep{
        time_step, std::move(reach_models),
        partial_evaluate(state_endpoint, n, initial_time)};
}

ValidatedMultiflowStep SelectedParametricMultifunctionIntegrator::flow_step(const ValidatedBoxMultifunctionInterface& field,
                                                                            const ExactBoxType& initial,
                                                                            const StepSizeType& time_step,
                                                                            const UpperBoxType& state_bound) const {
    return this->structured_step(field, initial, time_step, state_bound).bounding_step();
}

} // namespace Ariadne
