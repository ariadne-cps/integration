/***************************************************************************
 *            bounder.cpp
 *
 *  Copyright  2018-20  Luca Geretti
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

/** @file src/bounder.cpp
 *  @brief Implement bounder operations.
 */

#include "bounder.hpp"

#include "function/formula.hpp"
#include "function/taylor_model.hpp"

namespace Ariadne {
const MinimumStepSize DEFAULT_MINIMUM_STEP_SIZE = MinimumStepSize(0.00000095367431640625_x);
const Generator<MinimumStepSize> minimum_step_size = Generator<MinimumStepSize>();
const LipschitzTolerance DEFAULT_LIPSCHITZ_TOLERANCE = LipschitzTolerance(0.5_x);
const Generator<LipschitzTolerance> lipschitz_tolerance = Generator<LipschitzTolerance>();
Pair<StepSizeType, UpperBoxType> BounderHandle::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const Suggestion<StepSizeType>& hsug) const {
    return _impl->compute(f, D, hsug);
}
Pair<StepSizeType, UpperBoxType> BounderHandle::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const BoxDomainType& A, const Suggestion<StepSizeType>& hsug) const {
    return _impl->compute(f, D, A, hsug);
}
Pair<StepSizeType, UpperBoxType> BounderHandle::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const StepSizeType& t, const Suggestion<StepSizeType>& hsug) const {
    return _impl->compute(f, D, t, hsug);
}
Pair<StepSizeType, UpperBoxType> BounderHandle::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const StepSizeType& t, const BoxDomainType& A, const Suggestion<StepSizeType>& hsug) const {
    return _impl->compute(f, D, t, A, hsug);
}
BounderHandle* BounderHandle::clone() const { return new BounderHandle(*this); }
Void BounderHandle::_write(OutputStream& os) const { _impl->_write(os); }
BounderHandle::operator const BounderInterface&() const { return *_impl; }
BounderHandle& BounderHandle::operator=(const BounderHandle& other) {
    _impl = other._impl;
    return *this;
}
BounderHandle::BounderHandle(const BounderHandle& other) : _impl(other._impl) {}
BounderHandle::BounderHandle(const BounderInterface& bounder) : _impl(bounder.clone()) {}
EulerBounder* EulerBounder::clone() const { return new EulerBounder(*this); }
Void EulerBounder::_write(OutputStream& os) const { os << "EulerBounder"; }
MinimumStepSize::MinimumStepSize(const ExactDouble& value) : Attribute<ExactDouble>(value) {}
LipschitzTolerance::LipschitzTolerance(const ExactDouble& value) : Attribute<ExactDouble>(value) {}
BoundingNotFoundException::BoundingNotFoundException(const String& str) : std::runtime_error(str) {}

OutputStream& operator<<(OutputStream& os, const BounderInterface& bounder) {
    bounder._write(os);
    return os;
}

BounderBase::BounderBase(LipschitzTolerance lipschitz, MinimumStepSize minimum_step)
    : _lipschitz_tolerance(lipschitz), _minimum_step_size(minimum_step) {}

Pair<StepSizeType, UpperBoxType> BounderBase::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const Suggestion<StepSizeType>& hsug) const {
    return this->compute(f, D, BoxDomainType(0u), hsug);
}

Pair<StepSizeType, UpperBoxType> BounderBase::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const StepSizeType& t, const Suggestion<StepSizeType>& hsug) const {
    return this->compute(f, D, t, BoxDomainType(0u), hsug);
}

EulerBounder::EulerBounder(LipschitzTolerance lipschitz, MinimumStepSize minimum_step) : BounderBase(lipschitz, minimum_step) {}

Pair<StepSizeType, UpperBoxType> EulerBounder::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const BoxDomainType& A, const Suggestion<StepSizeType>& hsug) const {
    ARIADNE_PRECONDITION(f.result_size() == D.dimension());
    ARIADNE_PRECONDITION(f.argument_size() == D.dimension() + A.dimension());
    return this->_compute(f, D, 0, A, hsug);
}

Pair<StepSizeType, UpperBoxType> EulerBounder::compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const StepSizeType& t, const BoxDomainType& A, const Suggestion<StepSizeType>& hsug) const {
    ARIADNE_PRECONDITION(f.result_size() == D.dimension());
    ARIADNE_PRECONDITION(f.argument_size() == D.dimension() + 1u + A.dimension());
    return this->_compute(f, D, t, A, hsug);
}

Pair<StepSizeType, UpperBoxType> EulerBounder::_compute(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const StepSizeType& t, const BoxDomainType& A, const Suggestion<StepSizeType>& hsug) const {
    LOGGING_SCOPE_CREATE;
    const PositiveFloatDP BOX_RADIUS_WIDENING = cast_positive(0.25_exact);
    const PositiveFloatDP NO_WIDENING = cast_positive(1.0_exact);
    const PositiveFloatDP INITIAL_STARTING_WIDENING = cast_positive(2.0_exact);
    const PositiveFloatDP INITIAL_REFINING_WIDENING = cast_positive(1.125_exact);
    const StepSizeType MINIMUM_STEP_SIZE(1, 20u);
    const CounterType EXPANSION_STEPS = 4;
    const CounterType REFINEMENT_STEPS = 4;

    StepSizeType h = hsug.suggestion();

    FloatDPUpperBound lipschitz = norm(f.jacobian(Vector<FloatDPBounds>(cast_singleton(product(D, to_time_bounds(t, t + h), A))))).upper();
    StepSizeType hlip = static_cast<StepSizeType>(cast_exact(this->_lipschitz_tolerance.value() / lipschitz));
    h = min(hlip, h);
    LOGGING_PRINTLN("min(hlip,h)=" << h);

    IntervalDomainType T = to_time_bounds(t, t + h);

    LOGGING_PRINTLN("Finding contraction");

    UpperBoxType B = D;
    Bool success = false;
    StepSizeType hprev = h * 1.5_dy;
    while (!success) {
        B = this->_formula(f, D, T, A, D, BOX_RADIUS_WIDENING, INITIAL_STARTING_WIDENING);
        for (CounterType i = 0; i < EXPANSION_STEPS; ++i) {
            UpperBoxType Br = this->_refinement(f, D, T, A, B);
            if (not definitely(is_bounded(Br))) {
                success = false;
                LOGGING_PRINTLN_AT(1, "B is not bounded.");
                break;
            } else if (refines(Br, B)) {
                B = Br;
                success = true;
                LOGGING_PRINTLN_AT(1, "Found contraction of B=" << B);
                break;
            } else {
                B = this->_formula(f, D, T, A, B, NO_WIDENING, INITIAL_REFINING_WIDENING);
                LOGGING_PRINTLN_AT(1, "Expanding B to " << B);
            }
        }
        if (!success) {
            StepSizeType hnew = hlf(hprev);
            LOGGING_PRINTLN_AT(1, "Reduced h to " << h);
            hprev = h;
            h = StepSizeType(hnew.get_d());
            if (h < this->_minimum_step_size.value())
                ARIADNE_THROW(BoundingNotFoundException, "EulerBounder::_compute", "The step size is lower than the minimum (" << this->_minimum_step_size.value() << ") allowed, bounding could not be found.");
            T = to_time_bounds(t, t + h);
        }
    }

    LOGGING_PRINTLN("Refining B");
    for (CounterType i = 0; i < REFINEMENT_STEPS; ++i) {
        B = this->_refinement(f, D, T, A, B);
        LOGGING_PRINTLN_AT(1, "B=" << B);
    }

    LOGGING_PRINTLN("Found B=" << B << " using h=" << h);

    return std::make_pair(h, B);
}

UpperBoxType EulerBounder::_refinement(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const IntervalDomainType& T, const BoxDomainType& A, const UpperBoxType& B) const {
    const PositiveFloatDP NO_WIDENING = cast_positive(1.0_exact);
    return _formula(f, D, T, A, B, NO_WIDENING, NO_WIDENING);
}

UpperBoxType EulerBounder::_formula(const ValidatedVectorMultivariateFunction& f, const BoxDomainType& D, const IntervalDomainType& T, const BoxDomainType& A, const UpperBoxType& B, PositiveFloatDP INITIAL_BOX_WIDENING, PositiveFloatDP VECTOR_WIDENING) const {
    const UpperIntervalType rT(T);
    const UpperBoxType rA(A);
    const UpperIntervalType rH = rT - T.lower_bound();

    const bool is_autonomous = (f.argument_size() == D.dimension() + A.dimension());
    UpperBoxType dom = is_autonomous ? product(B, rA) : product(B, rT, rA);

    UpperBoxType wD = (INITIAL_BOX_WIDENING != 1) ? D : D + INITIAL_BOX_WIDENING * (D - D.midpoint());

    return wD + (VECTOR_WIDENING * rH) * cast_vector(apply(f, dom));
}

} // namespace Ariadne
