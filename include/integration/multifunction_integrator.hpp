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

#include <stdexcept>

#include "algebra/vector.hpp"
#include "function/taylor_multifunction.hpp"
#include "geometry/box.hpp"
#include "integration/integrator_interface.hpp"
#include "numeric/numeric.hpp"

namespace Ariadne {

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

/** @brief Box extension of a double-precision interval Taylor vector model. */
class ValidatedIntervalTaylorBoxMultifunctionDP : public ValidatedBoxMultifunctionInterface {
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

/** @brief First-order validated enclosure for x'(t) in F(x(t)).
 *
 *  Given B with X0 + [0,h] F(B) strictly inside B, continuity/continuation
 *  ensures every solution starting in X0 stays in B through time h.
 *  The method returns X0 + [0,h] F(B) for reach and X0 + h F(B) for evolve.
 *  The selection of uncertain coefficients is allowed to vary with time.
 *  No enclosure is returned when the strict continuation test fails.
 */
class FirstOrderMultifunctionIntegrator : public MultifunctionIntegratorInterface {
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
        if (not definitely(inside(reach, state_bound))) {
            throw FlowBoundsException("Cannot certify the multifunction reach tube inside the supplied state bound");
        }
        return ValidatedMultiflowStep{time_step, reach, evolve};
    }
};

} // namespace Ariadne

#endif /* ARIADNE_MULTIFUNCTION_INTEGRATOR_HPP */
