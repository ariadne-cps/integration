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
#include <variant>

#include "integration/integration_step_adapters.hpp"

namespace Ariadne {

/** @brief Validated step request, independent of the algorithm.
 *
 *  This is a numerical request, not a system model. The referenced field must
 *  outlive the call to step(). A continuation request retains the original
 *  state-coordinate domain and all earlier input parameters when a validated
 *  endpoint map is available. These numerical field overloads are autonomous.
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
    StepSizeType _initial_time;
    UpperBoxType _bounding_box;
    std::optional<ValidatedVectorMultivariateFunctionPatch> _initial_mapping;

  public:
    ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(initial), _time_step(time_step),
          _initial_time(0u), _bounding_box(bounding_box) {
        this->validate();
    }

    ValidatedIntegrationRequest(const ValidatedBoxMultifunctionInterface& field,
                                const ExactBoxType& initial,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(initial), _time_step(time_step),
          _initial_time(0u), _bounding_box(bounding_box) {
        this->validate();
    }

    /** @brief Derive the original state coordinates from a structured result,
     *  or conservatively restart with the full previous endpoint box.
     *  This does not select or retry a numerical algorithm.
     */
    static ExactBoxType continuation_domain(const ValidatedIntegrationStep& previous) {
        if (previous.has_evolve_model()) {
            const SizeType n = previous.evolve().dimension();
            if (previous.evolve_model().argument_size() < n) {
                throw std::invalid_argument("Previous endpoint mapping has too few state coordinates");
            }
            return ExactBoxType(previous.evolve_model().domain()[range(0u, n)]);
        }
        return cast_exact_box(previous.evolve());
    }

    /** @brief Continue from a validated step without inspecting its solver type.
     *  When previous has an endpoint map, retain it and its original domain.
     *  For a box-only result, the exact hull of the certified endpoint is used.
     */
    ValidatedIntegrationRequest(const ValidatedVectorMultivariateFunction& field,
                                const ValidatedIntegrationStep& previous,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(continuation_domain(previous)),
          _time_step(time_step), _initial_time(previous.final_time()),
          _bounding_box(bounding_box),
          _initial_mapping(previous.has_evolve_model()
                               ? std::optional<ValidatedVectorMultivariateFunctionPatch>(previous.evolve_model())
                               : std::nullopt) {
        this->validate();
    }

    ValidatedIntegrationRequest(const ValidatedBoxMultifunctionInterface& field,
                                const ValidatedIntegrationStep& previous,
                                const StepSizeType& time_step,
                                const UpperBoxType& bounding_box)
        : _field(std::cref(field)), _initial(continuation_domain(previous)),
          _time_step(time_step), _initial_time(previous.final_time()),
          _bounding_box(bounding_box),
          _initial_mapping(previous.has_evolve_model()
                               ? std::optional<ValidatedVectorMultivariateFunctionPatch>(previous.evolve_model())
                               : std::nullopt) {
        this->validate();
    }

    const Field& field() const { return _field; }
    const ExactBoxType& initial() const { return _initial; }
    const StepSizeType& time_step() const { return _time_step; }
    const StepSizeType& initial_time() const { return _initial_time; }
    const UpperBoxType& bounding_box() const { return _bounding_box; }
    bool has_initial_mapping() const { return _initial_mapping.has_value(); }
    const ValidatedVectorMultivariateFunctionPatch& initial_mapping() const {
        if (not _initial_mapping) {
            throw std::logic_error("This integration request has no initial mapping");
        }
        return *_initial_mapping;
    }

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
        if (_initial_mapping &&
            (_initial_mapping->result_size() != _initial.dimension() ||
             _initial_mapping->argument_size() < _initial.dimension() ||
             not (_initial_mapping->domain()[range(0u, _initial.dimension())] == _initial))) {
            throw std::invalid_argument("Validated continuation mapping has an incompatible initial-state domain");
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
        if (request.has_initial_mapping()) {
            throw std::invalid_argument("Ordinary flow continuation cannot discard a structured endpoint map");
        }
        const auto native = _integrator.flow_step(
            *field, request.initial(), request.time_step(), request.bounding_box());
        return validated_integration_step(native,
            Interval<StepSizeType>(StepSizeType(0u), request.time_step()),
            request.initial_time());
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
        if (request.has_initial_mapping()) {
            throw std::invalid_argument("Box-only integration cannot discard a structured endpoint map");
        }
        return validated_integration_step(_integrator.flow_step(
            *field, request.initial(), request.time_step(), request.bounding_box()),
            request.initial_time());
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
        if (request.has_initial_mapping()) {
            return validated_integration_step(_integrator.structured_step(
                *field, request.initial(), request.initial_mapping(),
                request.initial_time(), request.time_step(), request.bounding_box()),
                request.initial_time());
        }
        return validated_integration_step(_integrator.structured_step(
            *field, request.initial(), request.initial_time(),
            request.time_step(), request.bounding_box()),
            request.initial_time());
    }
};

} // namespace Ariadne

#endif /* ARIADNE_VALIDATED_STEP_INTEGRATOR_HPP */
