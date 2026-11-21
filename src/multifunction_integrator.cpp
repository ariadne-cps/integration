/***************************************************************************
 *            src/multifunction_integrator.cpp
 *
 *  Copyright  2026  Luca Geretti
 ****************************************************************************/

/** @file src/multifunction_integrator.cpp
 *  @brief Validated box multifunction adapter and first-order integration methods.
 */

#include "integration/multifunction_integrator.hpp"

namespace Ariadne {

ValidatedIntervalTaylorBoxMultifunctionDP::ValidatedIntervalTaylorBoxMultifunctionDP(const ValidatedVectorIntervalTaylorFunctionModel<FloatDP>& model)
    : _model(model) {}

SizeType ValidatedIntervalTaylorBoxMultifunctionDP::argument_size() const { return _model.argument_size(); }

SizeType ValidatedIntervalTaylorBoxMultifunctionDP::result_size() const { return _model.result_size(); }

BoxDomainType ValidatedIntervalTaylorBoxMultifunctionDP::domain() const { return _model.domain(); }

UpperBoxType ValidatedIntervalTaylorBoxMultifunctionDP::evaluate(const UpperBoxType& box) const {
    if (box.dimension() != this->argument_size() || not definitely(subset(box, this->domain()))) {
        throw std::invalid_argument("Multifunction box evaluation is outside its domain");
    }
    Vector<FloatDPBounds> arguments(box.dimension(), [&](SizeType i) {
        return FloatDPBounds(box[i].lower_bound(), box[i].upper_bound());
    });
    return _model(arguments);
}

InclusionIntegratorFieldData ValidatedIntervalTaylorBoxMultifunctionDP::inclusion_integrator_field() const {
    const auto polynomials = _model.polynomials();
    const SizeType n = _model.argument_size();
    List<IntervalDomainType> input_bounds;
    Vector<EffectiveFormula> expressions(_model.result_size(), EffectiveFormula::zero());

    for (SizeType i = 0u; i != _model.result_size(); ++i) {
        for (auto it = polynomials[i].begin(); it != polynomials[i].end(); ++it) {
            const auto& coefficient = it->coefficient();
            const SizeType parameter_index = n + input_bounds.size();
            input_bounds.append(IntervalDomainType(coefficient.lower_bound().raw(),
                                                   coefficient.upper_bound().raw()));

            EffectiveFormula monomial = EffectiveFormula::coordinate(parameter_index);
            for (SizeType j = 0u; j != n; ++j) {
                for (DegreeType p = 0u; p != it->index()[j]; ++p) {
                    monomial = monomial * EffectiveFormula::coordinate(j);
                }
            }
            expressions[i] = expressions[i] + monomial;
        }
    }
    if (input_bounds.empty()) {
        // Keep the structural representation nonempty even for the zero
        // polynomial, without adding any artificial uncertainty.
        input_bounds.append(IntervalDomainType(0.0_x, 0.0_x));
    }
    BoxDomainType inputs{Vector<IntervalDomainType>(input_bounds)};
    const SizeType arguments = n + inputs.dimension();
    return InclusionIntegratorFieldData{
        EffectiveVectorMultivariateFunction(
            EffectiveVectorMultivariateFormulaFunction(arguments, expressions)),
        inputs};
}

ValidatedTimedInputFieldBoxAdapterDP::ValidatedTimedInputFieldBoxAdapterDP(const EffectiveVectorMultivariateFunction& field,
                                                                           const BoxDomainType& states,
                                                                           const IntervalDomainType& time,
                                                                           const BoxDomainType& parameters,
                                                                           const BoxDomainType& inputs)
    : _field(field), _states(states), _time(time),
      _parameters(parameters), _inputs(inputs) {
    if (field.result_size() != states.dimension() ||
        field.argument_size() != states.dimension() + 1u +
                                     parameters.dimension() + inputs.dimension()) {
        throw std::invalid_argument("Timed multifunction field arguments must be [state,time,parameters,inputs]");
    }
}

SizeType ValidatedTimedInputFieldBoxAdapterDP::argument_size() const { return _states.dimension(); }

SizeType ValidatedTimedInputFieldBoxAdapterDP::result_size() const { return _field.result_size(); }

BoxDomainType ValidatedTimedInputFieldBoxAdapterDP::domain() const { return _states; }

IntervalDomainType ValidatedTimedInputFieldBoxAdapterDP::time_domain() const { return _time; }

BoxDomainType ValidatedTimedInputFieldBoxAdapterDP::parameter_domain() const { return _parameters; }

InclusionIntegratorFieldData ValidatedTimedInputFieldBoxAdapterDP::timed_inclusion_integrator_field() const {
    return InclusionIntegratorFieldData{_field, _inputs};
}

UpperBoxType ValidatedTimedInputFieldBoxAdapterDP::evaluate(const UpperBoxType& states, const UpperIntervalType& time,
                                                            const UpperBoxType& parameters) const {
    if (states.dimension() != _states.dimension() ||
        parameters.dimension() != _parameters.dimension() ||
        not definitely(subset(states, _states)) ||
        not definitely(subset(time, _time)) ||
        not definitely(subset(parameters, _parameters))) {
        throw std::invalid_argument("Timed multifunction evaluation outside its state, time or parameter domain");
    }
    const SizeType n = _states.dimension();
    const SizeType p = _parameters.dimension();
    Vector<FloatDPBounds> arguments(_field.argument_size(), [&](SizeType i) {
        if (i < n) {
            return FloatDPBounds(states[i].lower_bound(), states[i].upper_bound());
        }
        if (i == n) {
            return FloatDPBounds(time.lower_bound(), time.upper_bound());
        }
        if (i < n + 1u + p) {
            const auto& param = parameters[i - n - 1u];
            return FloatDPBounds(param.lower_bound(), param.upper_bound());
        }
        const auto& input = _inputs[i - n - 1u - p];
        return FloatDPBounds(input.lower_bound(), input.upper_bound());
    });
    const Vector<FloatDPBounds> values = _field(arguments);
    return UpperBoxType(values.size(), [&](SizeType i) {
        return UpperIntervalType(values[i].lower(), values[i].upper());
    });
}

ValidatedMultiflowStep FirstOrderMultifunctionIntegrator::flow_step(const ValidatedBoxMultifunctionInterface& field,
                                                                    const ExactBoxType& initial,
                                                                    const StepSizeType& time_step,
                                                                    const UpperBoxType& state_bound) const {
    if (not(time_step > StepSizeType(0u))) {
        throw std::invalid_argument("Multifunction integration requires a positive time step");
    }
    if (field.argument_size() != field.result_size() ||
        initial.dimension() != field.argument_size() ||
        state_bound.dimension() != field.argument_size()) {
        throw std::invalid_argument("Multifunction integration requires matching state and field dimensions");
    }
    if (not definitely(subset(state_bound, field.domain()))) {
        throw std::invalid_argument("Multifunction integration bound must lie inside the field domain");
    }

    const UpperBoxType values = field.evaluate(state_bound);
    if (values.dimension() != field.result_size()) {
        throw std::invalid_argument("Invalid multifunction box extension result dimension");
    }

    const UpperIntervalType duration(0u, time_step);
    const UpperIntervalType final_time(time_step);
    UpperBoxType reach(initial.dimension());
    UpperBoxType evolve(initial.dimension());
    for (SizeType i = 0; i != initial.dimension(); ++i) {
        const UpperIntervalType start(initial[i]);
        reach[i] = start + duration * values[i];
        evolve[i] = start + final_time * values[i];
    }
    // The supplied bounding box is defined by its exact floating-point
    // endpoints. Comparison between two outer-enclosure boxes alone may
    // be indeterminate even for strictly separated bounds.
    if (not definitely(inside(reach, cast_exact_box(state_bound)))) {
        std::ostringstream message;
        message << "Cannot certify multifunction tube: reach=" << reach << " bound=" << state_bound
                << " field_values=" << values << " step=" << time_step;
        throw FlowBoundsException(message.str());
    }
    return ValidatedMultiflowStep{time_step, reach, evolve};
}

ValidatedMultiflowStep FirstOrderMultifunctionIntegrator::flow_step(
    const ValidatedTimedBoxMultifunctionInterface& field,
    const ExactBoxType& initial,
    const StepSizeType& initial_time,
    const StepSizeType& time_step,
    const ExactBoxType& parameters,
    const UpperBoxType& state_bound) const {
    if (not(time_step > StepSizeType(0u))) {
        throw std::invalid_argument("Timed multifunction integration requires a positive time step");
    }
    if (field.argument_size() != field.result_size() ||
        initial.dimension() != field.argument_size() ||
        state_bound.dimension() != field.argument_size() ||
        parameters.dimension() != field.parameter_domain().dimension()) {
        throw std::invalid_argument("Timed multifunction integration has incompatible state or parameter dimensions");
    }
    const IntervalDomainType times(initial_time, initial_time + time_step);
    if (not definitely(subset(state_bound, field.domain())) ||
        not definitely(subset(times, field.time_domain())) ||
        not definitely(subset(parameters, field.parameter_domain()))) {
        throw std::invalid_argument("Timed multifunction step is outside its numerical field domain");
    }

    const UpperBoxType values = field.evaluate(
        state_bound, UpperIntervalType(initial_time, initial_time + time_step),
        UpperBoxType(parameters));
    if (values.dimension() != field.result_size()) {
        throw std::invalid_argument("Invalid timed multifunction box extension result dimension");
    }

    const UpperIntervalType duration(0u, time_step);
    const UpperIntervalType final_duration(time_step);
    UpperBoxType reach(initial.dimension());
    UpperBoxType evolve(initial.dimension());
    for (SizeType i = 0u; i != initial.dimension(); ++i) {
        const UpperIntervalType start(initial[i]);
        reach[i] = start + duration * values[i];
        evolve[i] = start + final_duration * values[i];
    }
    if (not definitely(inside(reach, cast_exact_box(state_bound)))) {
        std::ostringstream message;
        message << "Cannot certify timed multifunction tube: reach=" << reach
                << " bound=" << state_bound << " field_values=" << values
                << " time=[" << initial_time << "," << initial_time + time_step << "]";
        throw FlowBoundsException(message.str());
    }
    return ValidatedMultiflowStep{time_step, reach, evolve};
}

} // namespace Ariadne
