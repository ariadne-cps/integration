/***************************************************************************
 *            integration/multifunction_integrator.hpp
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

/** @file integration/multifunction_integrator.hpp
 *  @brief Validated first-order enclosures for differential inclusions.
 */

#ifndef ARIADNE_MULTIFUNCTION_INTEGRATOR_HPP
#define ARIADNE_MULTIFUNCTION_INTEGRATOR_HPP

#include <sstream>
#include <stdexcept>

#include "algebra/vector.hpp"
#include "function/function.hpp"
#include "function/symbolic_function.hpp"
#include "function/taylor_multifunction.hpp"
#include "geometry/box.hpp"
#include "integration/integrator_interface.hpp"
#include "numeric/numeric.hpp"

namespace Ariadne {


/** @brief Field data expected by the legacy high-order inclusion solver.
 *
 *  This is a numerical adapter payload, not a dynamical-system declaration.
 *  The field arguments are ordered as state coordinates followed by freely
 *  time-varying inputs. Trajectory-constant parameters are NOT supported here.
 *
 *  A conversion may conservatively enlarge the inclusion; independent interval
 *  coefficients may lose correlations, but must never remove trajectories.
 *  Availability of this payload does not select any integration algorithm.
 */
struct InclusionIntegratorFieldData {
    EffectiveVectorMultivariateFunction field;
    BoxDomainType inputs;
};

/** @brief A validated box extension of a set-valued vector field.
 *
 *  An implementation of evaluate(B) must enclose F(x) for every x in B.
 *  The vector field is autonomous, with compact convex nonempty values and
 *  sufficient regularity for local existence of inclusion solutions.
 *  This is a stronger contract than pointwise Multifunction evaluation.
 */
class ValidatedBoxMultifunctionInterface {
  public:
    virtual ~ValidatedBoxMultifunctionInterface() = default;
    virtual SizeType argument_size() const = 0;
    virtual SizeType result_size() const = 0;
    virtual BoxDomainType domain() const = 0;
    virtual UpperBoxType evaluate(const UpperBoxType& box) const = 0;

};

/** @brief Optional, solver-specific capability, independent of box evaluation.
 *
 *  This interface is implemented only by adapters that can produce certified
 *  f(x,u),U data for the legacy inclusion solver. It is not part of the base
 *  validated multifunction contract, does not declare a dynamical system and
 *  does not select an algorithm. Exposing the capability never freezes u(t).
 */
class InclusionIntegratorFieldProvider {
  public:
    virtual ~InclusionIntegratorFieldProvider() = default;
    virtual InclusionIntegratorFieldData inclusion_integrator_field() const = 0;
};

/** @brief Box extension of a double-precision interval Taylor vector model. */
class ValidatedIntervalTaylorBoxMultifunctionDP : public ValidatedBoxMultifunctionInterface,
                                                  public InclusionIntegratorFieldProvider {
    ValidatedVectorIntervalTaylorFunctionModel<FloatDP> _model;

  public:
    explicit ValidatedIntervalTaylorBoxMultifunctionDP(const ValidatedVectorIntervalTaylorFunctionModel<FloatDP>& model)
        : _model(model) {}

    virtual SizeType argument_size() const override { return _model.argument_size(); }
    virtual SizeType result_size() const override { return _model.result_size(); }
    virtual BoxDomainType domain() const override { return _model.domain(); }

    virtual UpperBoxType evaluate(const UpperBoxType& box) const override {
        if (box.dimension() != this->argument_size() || not definitely(subset(box, this->domain()))) {
            throw std::invalid_argument("Multifunction box evaluation is outside its domain");
        }
        Vector<FloatDPBounds> arguments(box.dimension(), [&](SizeType i) {
            return FloatDPBounds(box[i].lower_bound(), box[i].upper_bound());
        });
        return _model(arguments);
    }

    /** @brief Supply high-order solver data from an interval Taylor polynomial.
     *
     *  The public polynomials() method includes the Taylor remainder as an
     *  interval constant and applies the original-domain scaling. For every
     *  polynomial monomial c_alpha * x^alpha, a fresh input u_alpha in the
     *  (outward-rounded) coefficient interval is introduced. This gives an
     *  inclusion of the original pointwise multifunction, allowing u_alpha
     *  to vary over time, without treating it as a frozen ODE parameter.
     *
     *  This representation is conservative: distributing / unscaling the
     *  polynomial may widen bounds, and shared coefficient correlations
     *  are not reconstructed.
     */
    InclusionIntegratorFieldData inclusion_integrator_field() const override {
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
            inputs
        };
    }

};

/** @brief Rigorous reach-tube and endpoint box enclosures for a single step.
 *  Both components are certified by the integrator on successful return.
 *  The reach box contains all intermediate states, including the endpoints.
 */
struct ValidatedMultiflowStep {
    StepSizeType time_step;
    UpperBoxType reach;
    UpperBoxType evolve;
};

/** @brief Contract for validated integration of autonomous differential inclusions.
 *  Unlike IntegratorInterface, this returns independently meaningful reach
 *  and evolve enclosures, not a single ordinary-function flow patch.
 */
class MultifunctionIntegratorInterface {
  public:
    virtual ~MultifunctionIntegratorInterface() = default;
    virtual ValidatedMultiflowStep flow_step(const ValidatedBoxMultifunctionInterface& field,
                                             const ExactBoxType& initial,
                                             const StepSizeType& time_step,
                                             const UpperBoxType& state_bound) const = 0;
};


/** @brief Certified box extension for a time-dependent set-valued field.
 *
 *  The numerical field has state x, absolute time t and trajectory-constant
 *  parameters p, with arbitrary time-varying selections u(t) handled by
 *  evaluate(). For X, T and P, evaluate(X,T,P) must enclose every admissible
 *  instantaneous field value F(x,t,p), x in X, t in T, p in P.
 *
 *  This contract is weaker than a structure-preserving solver: it may lose
 *  dependence on p and time across steps, but may never exclude a trajectory.
 *  It does not introduce a dynamical-system declaration into integration.
 */
class ValidatedTimedBoxMultifunctionInterface {
  public:
    virtual ~ValidatedTimedBoxMultifunctionInterface() = default;
    virtual SizeType argument_size() const = 0;
    virtual SizeType result_size() const = 0;
    virtual BoxDomainType domain() const = 0;
    virtual IntervalDomainType time_domain() const = 0;
    virtual BoxDomainType parameter_domain() const = 0;
    virtual UpperBoxType evaluate(const UpperBoxType& states,
                                  const UpperIntervalType& time,
                                  const UpperBoxType& parameters) const = 0;
};

/** @brief Validated numerical adapter for f(x,t,p,u), with time-varying u.
 *
 *  Argument order is state, absolute time, fixed parameters, free inputs.
 *  The original function and domains are retained by value. Values of p stay
 *  constant along each admissible trajectory; interval box evaluation is
 *  intentionally conservative and does not encode their time correlation.
 *  Every u(t) in the input box remains admissible at each instant.
 *  This adapter does not supply high-order inclusion solver data.
 */
class ValidatedTimedInputFieldBoxAdapterDP : public ValidatedTimedBoxMultifunctionInterface {
    EffectiveVectorMultivariateFunction _field;
    BoxDomainType _states;
    IntervalDomainType _time;
    BoxDomainType _parameters;
    BoxDomainType _inputs;

  public:
    ValidatedTimedInputFieldBoxAdapterDP(const EffectiveVectorMultivariateFunction& field,
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

    SizeType argument_size() const override { return _states.dimension(); }
    SizeType result_size() const override { return _field.result_size(); }
    BoxDomainType domain() const override { return _states; }
    IntervalDomainType time_domain() const override { return _time; }
    BoxDomainType parameter_domain() const override { return _parameters; }

    UpperBoxType evaluate(const UpperBoxType& states, const UpperIntervalType& time,
                          const UpperBoxType& parameters) const override {
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
};

/** @brief Optional numerical solver capability for nonautonomous inclusions.
 *
 *  Integrators that only accept autonomous fields do not expose this method.
 *  The interface keeps each algorithm's applicability explicit.
 */
class TimedMultifunctionIntegratorInterface {
  public:
    virtual ~TimedMultifunctionIntegratorInterface() = default;
    virtual ValidatedMultiflowStep flow_step(
        const ValidatedTimedBoxMultifunctionInterface& field,
        const ExactBoxType& initial,
        const StepSizeType& initial_time,
        const StepSizeType& time_step,
        const ExactBoxType& parameters,
        const UpperBoxType& state_bound) const = 0;
};

/** @brief First-order validated enclosure for x'(t) in F(x(t)).
 *
 *  Given B with X0 + [0,h] F(B) strictly inside B, continuity/continuation
 *  ensures every solution starting in X0 stays in B through time h.
 *  The method returns X0 + [0,h] F(B) for reach and X0 + h F(B) for evolve.
 *  The selection of uncertain coefficients is allowed to vary with time.
 *  No enclosure is returned when the strict continuation test fails.
 */
class FirstOrderMultifunctionIntegrator : public MultifunctionIntegratorInterface,
                                           public TimedMultifunctionIntegratorInterface {
  public:
    virtual ValidatedMultiflowStep flow_step(const ValidatedBoxMultifunctionInterface& field,
                                             const ExactBoxType& initial,
                                             const StepSizeType& time_step,
                                             const UpperBoxType& state_bound) const override {
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

    /** @brief Validate an entire absolute-time step of x'(t) in F(x,t,p).
     *
     *  Evaluating over [t0,t0+h] and the entire fixed-parameter domain
     *  deliberately yields an outer enclosure. The strict tube-in-bound
     *  check gives a continuation certificate; it does not freeze u(t).
     */
    ValidatedMultiflowStep flow_step(
        const ValidatedTimedBoxMultifunctionInterface& field,
        const ExactBoxType& initial,
        const StepSizeType& initial_time,
        const StepSizeType& time_step,
        const ExactBoxType& parameters,
        const UpperBoxType& state_bound) const override {
        if (not (time_step > StepSizeType(0u))) {
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
};

} // namespace Ariadne

#endif /* ARIADNE_MULTIFUNCTION_INTEGRATOR_HPP */
