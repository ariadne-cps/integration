/***************************************************************************
 *            integration/selected_inclusion_integrator.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file integration/selected_inclusion_integrator.hpp
 *  @brief Explicit reach and evolve strategy selection for validated inclusions.
 */

#ifndef ARIADNE_SELECTED_INCLUSION_INTEGRATOR_HPP
#define ARIADNE_SELECTED_INCLUSION_INTEGRATOR_HPP

#include <stdexcept>
#include <utility>

#include "integration/inclusion_integrator.hpp"
#include "integration/integrator.hpp"
#include "integration/multifunction_integrator.hpp"

namespace Ariadne {

/** @brief Parametric high-order integration output with distinct reach and
 *  evolve representations.
 *
 *  The reach patches rigorously cover the complete trajectory tube.
 *  The evolve patch encloses the endpoint, possibly using another temporal
 *  input approximation. It must not be used as a reach-tube certificate.
 */
struct SelectedInclusionStep {
    StepSizeType time_step;
    List<ValidatedVectorMultivariateFunctionPatch> reach_models;
    ValidatedVectorMultivariateFunctionPatch evolve_model;

    /** @brief Discard patch structure to obtain the common validated box result.
     *  All reach patches are included in the hull, not selected by volume.
     *  The step duration is captured by the integrator, not relabelled by the caller.
     */
    ValidatedMultiflowStep bounding_step() const {
        if (reach_models.empty()) {
            throw FlowBoundsException("The selected reach method did not return a reach model");
        }
        UpperBoxType reach_box(reach_models[0u].range());
        for (SizeType i = 1u; i != reach_models.size(); ++i) {
            reach_box = hull(reach_box, UpperBoxType(reach_models[i].range()));
        }
        return ValidatedMultiflowStep{time_step, reach_box, UpperBoxType(evolve_model.range())};
    }
};

/** @brief Run exactly one configured reach method and one evolve method.
 *
 *  Zero and Constant are the current reach-rigorous input approximations.
 *  Any supported input approximation can be used for endpoint evolution.
 *  The selected methods are never ranked, retried, or replaced by another
 *  method when an integration step fails.
 *
 *  The parametric field and input domain are stored by value. The old
 *  InclusionIntegrator objects (which refer to their input data) are created
 *  only for the duration of step(), so this adapter is safely copyable.
 */
class SelectedInclusionIntegrator {
    EffectiveVectorMultivariateFunction _field;
    BoxDomainType _inputs;
    InclusionIntegratorFactory _factory;
    InputApproximation _reach_approximation;
    InputApproximation _evolve_approximation;

  public:
    SelectedInclusionIntegrator(const InclusionIntegratorFactory& factory,
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

    /** @brief Index of the one selected reach method. */
    Nat reach_index() const { return _reach_approximation.reference().index(); }
    /** @brief Index of the one selected endpoint method. */
    Nat evolve_index() const { return _evolve_approximation.reference().index(); }

    /** @brief Advance an existing validated state mapping over a time step.
     *
     *  The high-order reach is generated with the reach approximation.
     *  The endpoint is generated with the evolve approximation. When both
     *  selections agree, the reach model is reused without a second
     *  integration. Otherwise the evolve method computes its own model,
     *  which is used only at the end time.
     */
    SelectedInclusionStep step(const BoxDomainType& state_domain,
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
        if (not (initial_mapping.domain()[range(0u, n)] == state_domain)) {
            throw std::invalid_argument("Inclusion step initial mapping has an incompatible state domain");
        }
        auto reach_method = _factory.create(_field, _inputs, _reach_approximation);
        if (not reach_method.is_reach_rigorous()) {
            throw std::invalid_argument("The selected reach implementation does not certify the whole tube");
        }
        auto reaches = reach_method.reach(state_domain, initial_mapping, state_bound, initial_time, time_step);
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

        auto evolve_reaches = evolve_method.reach(state_domain, initial_mapping, state_bound, initial_time, time_step);
        if (evolve_reaches.size() != 1u) {
            throw FlowBoundsException("The selected evolve method did not return exactly one endpoint source");
        }
        auto endpoint = evolve_method.evolve(evolve_reaches[0u], final_time);
        return SelectedInclusionStep{time_step, std::move(reaches), std::move(endpoint)};
    }
};


/** @brief Numerical box-evaluation adapter for f(x,u), with input box U.
 *
 *  The structured function and U are retained by value. The interval
 *  evaluation encloses every allowed instantaneous value; u may change
 *  arbitrarily between times, rather than being fixed for the whole step.
 *  The supplied state domain bounds where this box extension may be used.
 */
class ValidatedInputFieldBoxAdapterDP : public ValidatedBoxMultifunctionInterface,
                                        public InclusionIntegratorFieldProvider {
    EffectiveVectorMultivariateFunction _field;
    BoxDomainType _inputs;
    BoxDomainType _state_domain;

  public:
    ValidatedInputFieldBoxAdapterDP(const EffectiveVectorMultivariateFunction& field,
                                         const BoxDomainType& inputs,
                                         const BoxDomainType& state_domain)
        : _field(field), _inputs(inputs), _state_domain(state_domain) {
        if (inputs.dimension() == 0u || field.result_size() != state_domain.dimension() ||
            field.argument_size() != field.result_size() + inputs.dimension()) {
            throw std::invalid_argument("Parametric multifunction requires n states and at least one input");
        }
    }

    SizeType argument_size() const override { return _field.result_size(); }
    SizeType result_size() const override { return _field.result_size(); }
    BoxDomainType domain() const override { return _state_domain; }

    /** @brief Structural ODE-with-inputs representation for high-order solvers. */
    const EffectiveVectorMultivariateFunction& vector_field() const { return _field; }
    /** @brief Independent admissible input box, without freezing inputs in time. */
    const BoxDomainType& inputs() const { return _inputs; }

    InclusionIntegratorFieldData inclusion_integrator_field() const override {
        return InclusionIntegratorFieldData{_field, _inputs};
    }

    UpperBoxType evaluate(const UpperBoxType& box) const override {
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
};

/** @brief High-order inclusion integrator implementing the common box interface.
 *
 *  Higher-order integration requires explicitly supplied solver field data.
 *  Either a structured f(x,u) input or a Taylor model with interval coefficients
 *  can provide this data by implementing a separate numerical capability.
 *  A solver must reject fields that do not provide this capability rather than
 *  inventing unsupported correlations.
 *
 *  The optional structured step retains parameter correlations for downstream
 *  evolution and reconditioning. flow_step() exposes the common box result.
 */
class SelectedParametricMultifunctionIntegrator : public MultifunctionIntegratorInterface {
    InclusionIntegratorFactory _factory;
    ValidatedFunctionPatchFactory _function_factory;
    InputApproximation _reach_approximation;
    InputApproximation _evolve_approximation;

  public:
    SelectedParametricMultifunctionIntegrator(const IntegratorBase& integrator,
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

    Nat reach_index() const { return _reach_approximation.reference().index(); }
    Nat evolve_index() const { return _evolve_approximation.reference().index(); }

    /** @brief Start a validated structured inclusion step at time zero. */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
                                          const ExactBoxType& initial,
                                          const StepSizeType& time_step,
                                          const UpperBoxType& state_bound) const {
        return this->structured_step(field, initial, TimeStepType(0u),
                                     time_step, state_bound);
    }

    /** @brief Start from a fresh initial-state identity at an absolute time. */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
                                          const ExactBoxType& initial,
                                          const TimeStepType& initial_time,
                                          const StepSizeType& time_step,
                                          const UpperBoxType& state_bound) const {
        const auto identity = _function_factory.create_identity(initial);
        return this->structured_step(field, initial, identity,
                                     initial_time, time_step, state_bound);
    }

    /** @brief Continue a prior validated endpoint map without resetting
     *  its original state coordinates or previously introduced input parameters.
     */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
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
            not (initial_mapping.domain()[range(0u, initial.dimension())] == initial)) {
            throw std::invalid_argument("Structured continuation has incompatible mapping coordinates");
        }
        if (not definitely(subset(UpperBoxType(initial_mapping.range()),
                                  cast_exact_box(state_bound)))) {
            throw std::invalid_argument("Structured continuation endpoint is outside the certified state bound");
        }
        SelectedInclusionIntegrator selected(_factory, form.field, form.inputs,
                                            _reach_approximation, _evolve_approximation);
        return selected.step(initial, initial_mapping, state_bound, initial_time, time_step);
    }


    /** @brief Rigorous structured integration of x' in F(x,t,p,u(t)).
     *
     *  A numerical autonomous lifting defines (x,clock,p)'=(f(x,clock,p,u),1,0).
     *  Thus fixed p and absolute time are integrated as STATES, never as
     *  freely varying solver inputs. The only free inputs remain u(t).
     *  Output maps retain x, p and old input parameters; the synthetic
     *  initial-clock coordinate is removed by validated partial evaluation.
     */
    SelectedInclusionStep timed_structured_step(
        const ValidatedTimedBoxMultifunctionInterface& field,
        const ExactBoxType& initial,
        const TimeStepType& initial_time,
        const StepSizeType& time_step,
        const ExactBoxType& parameters,
        const UpperBoxType& state_bound,
        const ValidatedVectorMultivariateFunctionPatch* previous = nullptr) const {
        const SizeType n = initial.dimension();
        const SizeType np = parameters.dimension();
        const SizeType ns = n + 1u + np;
        const IntervalDomainType times(initial_time, initial_time + time_step);
        if (not (time_step > StepSizeType(0u)) || field.argument_size() != n ||
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
            return i < n ? original[i] :
                (i == n ? EffectiveFormula::constant(EffectiveNumber(1u))
                        : EffectiveFormula::zero());
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
                not (previous->domain()[range(0u, n)] == BoxDomainType(initial)) ||
                not (previous->domain()[range(n, n + np)] == BoxDomainType(parameters))) {
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

    /** @brief Compute certified boxes through the common integrator contract. */
    ValidatedMultiflowStep flow_step(const ValidatedBoxMultifunctionInterface& field,
                                    const ExactBoxType& initial,
                                    const StepSizeType& time_step,
                                    const UpperBoxType& state_bound) const override {
        return this->structured_step(field, initial, time_step, state_bound).bounding_step();
    }
};

} // namespace Ariadne

#endif // ARIADNE_SELECTED_INCLUSION_INTEGRATOR_HPP
