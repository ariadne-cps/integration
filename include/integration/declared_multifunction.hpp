/***************************************************************************
 *            integration/declared_multifunction.hpp
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

/** @file integration/declared_multifunction.hpp
 *  @brief Validated box field with separately declared constant and time-varying inputs.
 */

#ifndef ARIADNE_DECLARED_MULTIFUNCTION_HPP
#define ARIADNE_DECLARED_MULTIFUNCTION_HPP

#include <optional>
#include <stdexcept>

#include "integration/multifunction_integrator.hpp"

namespace Ariadne {

/** @brief Declarative autonomous dynamics x'=f(x,p,u(t)).
 *
 *  The coordinates of the formula are ordered as state x, trajectory-constant
 *  parameters p, then freely time-varying inputs u. A parameter p is selected
 *  once for each trajectory and remains unchanged during all subsequent steps.
 *  The box evaluation deliberately forgets this temporal dependency and is a
 *  sound, possibly larger, instantaneous enclosure.
 *
 *  No integrator is selected by this representation.
 */
class ValidatedDeclaredBoxMultifunctionDP : public ValidatedBoxMultifunctionInterface {
    EffectiveVectorMultivariateFunction _field;
    BoxDomainType _state_domain;
    BoxDomainType _constant_parameters;
    BoxDomainType _time_inputs;

  public:
    ValidatedDeclaredBoxMultifunctionDP(const EffectiveVectorMultivariateFunction& field,
                                        const BoxDomainType& state_domain,
                                        const BoxDomainType& constant_parameters,
                                        const BoxDomainType& time_inputs)
        : _field(field),
          _state_domain(state_domain),
          _constant_parameters(constant_parameters),
          _time_inputs(time_inputs) {
        if (state_domain.dimension() == 0u ||
            field.result_size() != state_domain.dimension() ||
            field.argument_size() != state_domain.dimension() +
                                     constant_parameters.dimension() +
                                     time_inputs.dimension()) {
            throw std::invalid_argument(
                "Declared multifunction requires state, parameter and input dimensions matching the field");
        }
    }

    SizeType argument_size() const override { return _state_domain.dimension(); }
    SizeType result_size() const override { return _field.result_size(); }
    BoxDomainType domain() const override { return _state_domain; }

    const EffectiveVectorMultivariateFunction& vector_field() const { return _field; }
    const BoxDomainType& constant_parameters() const { return _constant_parameters; }
    const BoxDomainType& time_inputs() const { return _time_inputs; }

    /** @brief Rigorous instantaneous range on a state box.
     *
     *  Encloses f(x,p,u) for every x in box, p in the constant-parameter
     *  domain and u in the time-varying-input domain. It does not assert
     *  that p can change along a single trajectory.
     */
    UpperBoxType evaluate(const UpperBoxType& box) const override {
        if (box.dimension() != this->argument_size() ||
            not definitely(subset(box, _state_domain))) {
            throw std::invalid_argument("Declared multifunction evaluation outside the state domain");
        }
        const SizeType n = _state_domain.dimension();
        const SizeType np = _constant_parameters.dimension();
        const SizeType nu = _time_inputs.dimension();
        Vector<FloatDPBounds> arguments(n + np + nu, [&](SizeType i) {
            if (i < n) {
                return FloatDPBounds(box[i].lower_bound(), box[i].upper_bound());
            }
            if (i < n + np) {
                const auto& parameter = _constant_parameters[i - n];
                return FloatDPBounds(parameter.lower_bound(), parameter.upper_bound());
            }
            const auto& input = _time_inputs[i - n - np];
            return FloatDPBounds(input.lower_bound(), input.upper_bound());
        });
        const Vector<FloatDPBounds> values = _field(arguments);
        return UpperBoxType(values.size(), [&](SizeType i) {
            return UpperIntervalType(values[i].lower(), values[i].upper());
        });
    }

    /** @brief Expose the ordinary inclusion backend only if no fixed
     *  parameter can accidentally be treated as a time-varying input.
     *
     *  A parameter-free ODE has no input coordinates and also does not
     *  supply this capability: it belongs to the ordinary ODE backend.
     */
    std::optional<ParametricMultifunctionForm> parametric_form() const override {
        if (_constant_parameters.dimension() != 0u || _time_inputs.dimension() == 0u) {
            return std::nullopt;
        }
        return ParametricMultifunctionForm{_field, _time_inputs};
    }

    /** @brief Explicitly relax trajectory-constant parameters into
     *  freely time-varying inputs, retaining their coordinate sharing.
     *
     *  The resulting differential inclusion contains every trajectory
     *  of the original parameterised system, but usually admits more.
     *  Callers must explicitly opt into this loss of temporal dependency.
     */
    ParametricMultifunctionForm time_varying_relaxation() const {
        if (_constant_parameters.dimension() + _time_inputs.dimension() == 0u) {
            throw std::invalid_argument("A parameter-free ODE has no inputs to relax");
        }
        List<IntervalDomainType> bounds;
        for (SizeType i = 0u; i != _constant_parameters.dimension(); ++i) {
            bounds.append(_constant_parameters[i]);
        }
        for (SizeType i = 0u; i != _time_inputs.dimension(); ++i) {
            bounds.append(_time_inputs[i]);
        }
        BoxDomainType relaxed_inputs{Vector<IntervalDomainType>(bounds)};
        return ParametricMultifunctionForm{_field, relaxed_inputs};
    }
};

} // namespace Ariadne

#endif /* ARIADNE_DECLARED_MULTIFUNCTION_HPP */
