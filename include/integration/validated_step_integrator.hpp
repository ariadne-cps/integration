/***************************************************************************
 *            integration/validated_step_integrator.hpp
 *
 *  Copyright  2026  Luca Geretti
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file integration/validated_step_integrator.hpp
 *  @brief Common numerical step interface for different validated field types.
 */

#ifndef ARIADNE_VALIDATED_STEP_INTEGRATOR_HPP
#define ARIADNE_VALIDATED_STEP_INTEGRATOR_HPP

#include <functional>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <variant>

#include "integration/integration_step_adapters.hpp"

namespace Ariadne {

/** @brief Validated step request, independent of the algorithm.
 *
 *  This is a numerical request, not a system model. Single-valued function
 *  handles are retained by value, including implicit conversions from effective
 *  functions. Polymorphic multifunction fields are borrowed from lvalues and
 *  must outlive every request referring to them. Temporaries cannot be borrowed.
 *  A continuation request retains the original
 *  state-coordinate domain and all earlier input parameters when a validated
 *  endpoint map is available. Numerical field overloads support both
 *  autonomous and explicit-time requests, according to solver capability.
 *
 *  A single-valued function can still represent a system with uncertain
 *  parameters: the distinction here is the numerical input capability,
 *  NOT the classification of the dynamical system.
 */
class ValidatedIntegrationRequest {
  public:
    using Field = std::variant<
        ValidatedVectorMultivariateFunction,
        std::reference_wrapper<const ValidatedBoxMultifunctionInterface>,
        std::reference_wrapper<const ValidatedTimedBoxMultifunctionInterface>>;

  private:
    Field _field;
    ExactBoxType _initial;
    StepSizeType _time_step;
    StepSizeType _initial_time;
    UpperBoxType _bounding_box;
    std::optional<ValidatedVectorMultivariateFunctionPatch> _initial_mapping;
    bool _has_explicit_time = false;
    ExactBoxType _constant_parameter_domain = ExactBoxType(0u);

  public:
ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box);


ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedBoxMultifunctionInterface> field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box);


    /** @brief Fresh nonautonomous equation step with fixed parameters.
     *  The field arguments are [state, absolute time, constant parameters].
     *  A zero-dimensional parameter domain is permitted.
     */
ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ExactBoxType& initial,
                                const StepSizeType& initial_time,
                                const StepSizeType& time_step,
                                const ExactBoxType& constant_parameters,
                                const UpperBoxType& bounding_box);


    /** @brief Fresh timed inclusion step with fixed parameter coordinates.
     *  Input uncertainty is represented by the numerical field itself.
     */
ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedTimedBoxMultifunctionInterface> field,
                                const ExactBoxType& initial,
                                const StepSizeType& initial_time,
                                const StepSizeType& time_step,
                                const ExactBoxType& constant_parameters,
                                const UpperBoxType& bounding_box);


    /** @brief Derive the original state coordinates from a structured result,
     *  or conservatively restart with the full previous endpoint box.
     *  This does not select or retry a numerical algorithm.
     */
    static ExactBoxType continuation_domain(const ValidatedIntegrationStep& previous);

    /** @brief Continue from a validated step without inspecting its solver type.
     *  When previous has an endpoint map, retain it and its original domain.
     *  For a box-only result, the exact hull of the certified endpoint is used.
     */
ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ValidatedIntegrationStep& previous,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box);


ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedBoxMultifunctionInterface> field,
                                const ValidatedIntegrationStep& previous,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box);


    /** @brief Continue a timed inclusion from the previous common result.
     *  Its parameter domain is retained, not independently reselected.
     *  A box-only algorithm still cannot accept a structured endpoint map.
     */
ValidatedIntegrationRequest(std::reference_wrapper<const ValidatedTimedBoxMultifunctionInterface> field,
                                const ValidatedIntegrationStep& previous,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box);


    const Field& field() const;
    const ExactBoxType& initial() const;
    const StepSizeType& time_step() const;
    const StepSizeType& initial_time() const;
    const UpperBoxType& bounding_box() const;
    /** @brief Generic numerical metadata, also valid for timed inclusions. */
    bool has_explicit_time() const;
    const ExactBoxType& constant_parameter_domain() const;
    bool has_initial_mapping() const;
    const ValidatedVectorMultivariateFunctionPatch& initial_mapping() const;

    template<class T>
    const T* field_as() const {
        if constexpr (std::is_same_v<T, ValidatedVectorMultivariateFunction>) {
            return std::get_if<T>(&_field);
        } else {
            const auto* reference = std::get_if<std::reference_wrapper<const T>>(&_field);
            return reference == nullptr ? nullptr : &reference->get();
        }
    }

  private:
    void validate() const;
};

/** @brief Common step operation for all supported validated solvers.
 *
 *  A concrete adapter accepts only capabilities the underlying algorithm
 *  implements. Unsupported numerical inputs fail explicitly; no algorithm
 *  is selected, retried or substituted by this interface.
 */
class ValidatedStepIntegratorInterface {
  public:
    virtual ~ValidatedStepIntegratorInterface() = default;
    virtual ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const = 0;
};

/** @brief Bridge the existing ordinary flow API to the common result. */
class FlowStepIntegratorAdapter : public ValidatedStepIntegratorInterface {
    const IntegratorInterface& _integrator;

  public:
    explicit FlowStepIntegratorAdapter(const IntegratorInterface& integrator);


    /** @brief The borrowed solver must outlive this adapter; temporaries are rejected. */
    FlowStepIntegratorAdapter(const IntegratorInterface&&) = delete;

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override;
};

/** @brief Bridge an inclusion box integrator to the common result. */
class MultiflowStepIntegratorAdapter : public ValidatedStepIntegratorInterface {
    const MultifunctionIntegratorInterface& _integrator;

  public:
    explicit MultiflowStepIntegratorAdapter(const MultifunctionIntegratorInterface& integrator);


    /** @brief The borrowed solver must outlive this adapter; temporaries are rejected. */
    MultiflowStepIntegratorAdapter(const MultifunctionIntegratorInterface&&) = delete;

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override;
};

/** @brief Retain high-order reach and endpoint patches via the common result. */
class StructuredInclusionStepIntegratorAdapter : public ValidatedStepIntegratorInterface {
    const SelectedParametricMultifunctionIntegrator& _integrator;

  public:
    explicit StructuredInclusionStepIntegratorAdapter(
        const SelectedParametricMultifunctionIntegrator& integrator);


    /** @brief The borrowed solver must outlive this adapter; temporaries are rejected. */
    StructuredInclusionStepIntegratorAdapter(const SelectedParametricMultifunctionIntegrator&&) = delete;

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override;
};

} // namespace Ariadne

#endif /* ARIADNE_VALIDATED_STEP_INTEGRATOR_HPP */
