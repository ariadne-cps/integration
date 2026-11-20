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
    ValidatedMultiflowStep bounding_step() const;
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
                                const InputApproximation& evolve_approximation);


    /** @brief Index of the one selected reach method. */
    Nat reach_index() const;
    /** @brief Index of the one selected endpoint method. */
    Nat evolve_index() const;

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
                               const StepSizeType& time_step) const;
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
                                         const BoxDomainType& state_domain);


    SizeType argument_size() const override;
    SizeType result_size() const override;
    BoxDomainType domain() const override;

    /** @brief Structural ODE-with-inputs representation for high-order solvers. */
    const EffectiveVectorMultivariateFunction& vector_field() const;
    /** @brief Independent admissible input box, without freezing inputs in time. */
    const BoxDomainType& inputs() const;

    InclusionIntegratorFieldData inclusion_integrator_field() const override;

    UpperBoxType evaluate(const UpperBoxType& box) const override;
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
                                             const InputApproximation& evolve_approximation);


    Nat reach_index() const;
    Nat evolve_index() const;

    /** @brief Start a validated structured inclusion step at time zero. */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
                                          const ExactBoxType& initial,
                                          const StepSizeType& time_step,
                                          const UpperBoxType& state_bound) const;

    /** @brief Start from a fresh initial-state identity at an absolute time. */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
                                          const ExactBoxType& initial,
                                          const TimeStepType& initial_time,
                                          const StepSizeType& time_step,
                                          const UpperBoxType& state_bound) const;

    /** @brief Continue a prior validated endpoint map without resetting
     *  its original state coordinates or previously introduced input parameters.
     */
    SelectedInclusionStep structured_step(const ValidatedBoxMultifunctionInterface& field,
                                          const ExactBoxType& initial,
                                          const ValidatedVectorMultivariateFunctionPatch& initial_mapping,
                                          const TimeStepType& initial_time,
                                          const StepSizeType& time_step,
                                          const UpperBoxType& state_bound) const;


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
        const ValidatedVectorMultivariateFunctionPatch* previous = nullptr) const;

    /** @brief Compute certified boxes through the common integrator contract. */
    ValidatedMultiflowStep flow_step(const ValidatedBoxMultifunctionInterface& field,
                                    const ExactBoxType& initial,
                                    const StepSizeType& time_step,
                                    const UpperBoxType& state_bound) const override;
};

} // namespace Ariadne

#endif // ARIADNE_SELECTED_INCLUSION_INTEGRATOR_HPP
