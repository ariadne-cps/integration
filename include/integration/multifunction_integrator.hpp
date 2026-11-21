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

/** @brief Optional structured numerical capability for nonautonomous inclusions.
 *  The field arguments are [state, absolute time, fixed parameters, free inputs].
 *  Only the final block contains independently selectable functions u(t).
 */
class TimedInclusionIntegratorFieldProvider {
  public:
    virtual ~TimedInclusionIntegratorFieldProvider() = default;
    virtual InclusionIntegratorFieldData timed_inclusion_integrator_field() const = 0;
};

/** @brief Box extension of a double-precision interval Taylor vector model. */
class ValidatedIntervalTaylorBoxMultifunctionDP : public ValidatedBoxMultifunctionInterface,
                                                  public InclusionIntegratorFieldProvider {
    ValidatedVectorIntervalTaylorFunctionModel<FloatDP> _model;

  public:
    explicit ValidatedIntervalTaylorBoxMultifunctionDP(const ValidatedVectorIntervalTaylorFunctionModel<FloatDP>& model);

    virtual SizeType argument_size() const override;
    virtual SizeType result_size() const override;
    virtual BoxDomainType domain() const override;

    virtual UpperBoxType evaluate(const UpperBoxType& box) const override;

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
    InclusionIntegratorFieldData inclusion_integrator_field() const override;
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
class ValidatedTimedInputFieldBoxAdapterDP : public ValidatedTimedBoxMultifunctionInterface,
                                             public TimedInclusionIntegratorFieldProvider {
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
                                         const BoxDomainType& inputs);

    SizeType argument_size() const override;
    SizeType result_size() const override;
    BoxDomainType domain() const override;
    IntervalDomainType time_domain() const override;
    BoxDomainType parameter_domain() const override;

    InclusionIntegratorFieldData timed_inclusion_integrator_field() const override;

    UpperBoxType evaluate(const UpperBoxType& states, const UpperIntervalType& time,
                          const UpperBoxType& parameters) const override;
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
                                             const UpperBoxType& state_bound) const override;

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
        const UpperBoxType& state_bound) const override;
};

} // namespace Ariadne

#endif /* ARIADNE_MULTIFUNCTION_INTEGRATOR_HPP */
