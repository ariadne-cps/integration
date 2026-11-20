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
    List<ValidatedVectorMultivariateFunctionPatch> reach_models;
    ValidatedVectorMultivariateFunctionPatch evolve_model;

    /** @brief Discard patch structure to obtain the common validated box result.
     *  All reach patches are included in the hull, not selected by volume.
     */
    ValidatedMultiflowStep bounding_step(const StepSizeType& time_step) const {
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
            return SelectedInclusionStep{std::move(reaches), std::move(endpoint)};
        }

        auto evolve_reaches = evolve_method.reach(state_domain, initial_mapping, state_bound, initial_time, time_step);
        if (evolve_reaches.size() != 1u) {
            throw FlowBoundsException("The selected evolve method did not return exactly one endpoint source");
        }
        auto endpoint = evolve_method.evolve(evolve_reaches[0u], final_time);
        return SelectedInclusionStep{std::move(reaches), std::move(endpoint)};
    }
};

} // namespace Ariadne

#endif // ARIADNE_SELECTED_INCLUSION_INTEGRATOR_HPP
