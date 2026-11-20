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
#include <stdexcept>
#include <variant>

#include "integration/integration_step_adapters.hpp"

namespace Ariadne {

/** @brief Autonomous validated step request, independent of the algorithm.
 *
 *  This is a numerical request, not a system model. The referenced field and
 *  initial/bounding boxes must outlive the call to step(). In this first bridge
 *  the time origin is zero; nonautonomous fields and multi-step integration
 *  remain available through the existing specialized interfaces.
 *
 *  A single-valued function can still represent a system with uncertain
 *  parameters: the distinction here is the numerical input capability,
 *  NOT the classification of the dynamical system.
 */
class ValidatedIntegrationRequest {
  public:
    using Field = std::variant<
        std::reference_wrapper<const ValidatedVectorMultivariateFunction>,
        std::reference_wrapper<const ValidatedBoxMultifunctionInterface>>;

  private:
    Field _field;
    ExactBoxType _initial;
    StepSizeType _time_step;
    UpperBoxType _bounding_box;

  public:
    ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(initial), _time_step(time_step),
          _bounding_box(bounding_box) {
        this->validate();
    }

    ValidatedIntegrationRequest(const ValidatedBoxMultifunctionInterface& field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(initial), _time_step(time_step),
          _bounding_box(bounding_box) {
        this->validate();
    }

    const Field& field() const { return _field; }
    const ExactBoxType& initial() const { return _initial; }
    const StepSizeType& time_step() const { return _time_step; }
    const UpperBoxType& bounding_box() const { return _bounding_box; }

    template<class T>
    const T* field_as() const {
        const auto* reference = std::get_if<std::reference_wrapper<const T>>(&_field);
        return reference == nullptr ? nullptr : &reference->get();
    }

  private:
    void validate() const {
        if (not (_time_step > StepSizeType(0u)) ||
            _initial.dimension() != _bounding_box.dimension()) {
            throw std::invalid_argument("Validated integration request has incompatible dimensions or duration");
        }
    }
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
    explicit FlowStepIntegratorAdapter(const IntegratorInterface& integrator)
        : _integrator(integrator) {}

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override {
        const auto* field = request.field_as<ValidatedVectorMultivariateFunction>();
        if (field == nullptr) {
            throw std::invalid_argument("This flow integrator requires a validated single-valued vector field");
        }
        const auto native = _integrator.flow_step(
            *field, request.initial(), request.time_step(), request.bounding_box());
        return validated_integration_step(native,
            Interval<StepSizeType>(StepSizeType(0u), request.time_step()));
    }
};

/** @brief Bridge an inclusion box integrator to the common result. */
class MultiflowStepIntegratorAdapter : public ValidatedStepIntegratorInterface {
    const MultifunctionIntegratorInterface& _integrator;

  public:
    explicit MultiflowStepIntegratorAdapter(const MultifunctionIntegratorInterface& integrator)
        : _integrator(integrator) {}

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override {
        const auto* field = request.field_as<ValidatedBoxMultifunctionInterface>();
        if (field == nullptr) {
            throw std::invalid_argument("This inclusion integrator requires a validated box multifunction");
        }
        return validated_integration_step(_integrator.flow_step(
            *field, request.initial(), request.time_step(), request.bounding_box()));
    }
};

/** @brief Retain high-order reach and endpoint patches via the common result. */
class StructuredInclusionStepIntegratorAdapter : public ValidatedStepIntegratorInterface {
    const SelectedParametricMultifunctionIntegrator& _integrator;

  public:
    explicit StructuredInclusionStepIntegratorAdapter(
        const SelectedParametricMultifunctionIntegrator& integrator)
        : _integrator(integrator) {}

    ValidatedIntegrationStep step(const ValidatedIntegrationRequest& request) const override {
        const auto* field = request.field_as<ValidatedBoxMultifunctionInterface>();
        if (field == nullptr) {
            throw std::invalid_argument("Structured inclusion integration requires a validated box multifunction");
        }
        return validated_integration_step(_integrator.structured_step(
            *field, request.initial(), request.time_step(), request.bounding_box()));
    }
};

} // namespace Ariadne

#endif /* ARIADNE_VALIDATED_STEP_INTEGRATOR_HPP */
