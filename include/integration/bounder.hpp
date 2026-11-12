/***************************************************************************
 *            integration/bounder.hpp
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

/** @file integration/bounder.hpp
 *  @brief Flow bounding strategies and interfaces for ordinary differential equations.
 */

#ifndef ARIADNE_BOUNDER_HPP
#define ARIADNE_BOUNDER_HPP

#include "utility/typedefs.hpp"
#include "utility/attribute.hpp"
#include "algebra/algebra.hpp"
#include "function/domain.hpp"
#include "function/function_model.hpp"
#include "logging/logging.hpp"
#include "integrator_interface.hpp"


namespace Ariadne {


//! \brief An argument is a suggestion or hint, which is allowed to be ignored or used approximately.
template<class T> class Suggestion {
    T _suggestion;
  public:
    explicit Suggestion(T const& t) : _suggestion(t) { }
    template<class TT> requires Convertible<TT,T> Suggestion(Suggestion<TT> const& s)
        : Suggestion(s.suggestion()) { }
    /** @brief Perform the integration operation. */
    template<class TT> friend Suggestion<TT> suggest(TT&& t);
    T const& suggestion() const { return _suggestion; }
    explicit operator T& () { return _suggestion; }
    explicit operator T const& () const { return _suggestion; }
};
template<class T> Suggestion<T> suggest(T&& t) { return Suggestion<T>{std::forward<T>(t)}; }

/** @brief Error raised when a flow bound cannot be established. */
class BoundingNotFoundException : public std::runtime_error {
  public:
    BoundingNotFoundException(const String& str) : std::runtime_error(str) { }
};

/** @brief Type representing lipschitz tolerance in integration algorithms. */
struct LipschitzTolerance : Attribute<ExactDouble> {
    explicit LipschitzTolerance(ExactDouble const& value) : Attribute<ExactDouble>(value) { }
};
static const Generator<LipschitzTolerance> lipschitz_tolerance = Generator<LipschitzTolerance>();
/** @brief Construct the integration object. */
static const LipschitzTolerance DEFAULT_LIPSCHITZ_TOLERANCE(0.5_x);
/** @brief Type representing minimum step size in integration algorithms. */
struct MinimumStepSize : Attribute<ExactDouble> {
    explicit MinimumStepSize(ExactDouble const& value) : Attribute<ExactDouble>(value) { }
};
static const Generator<MinimumStepSize> minimum_step_size = Generator<MinimumStepSize>();
/** @brief Construct the integration object. */
static const MinimumStepSize DEFAULT_MINIMUM_STEP_SIZE(0.00000095367431640625_x);

//! \ingroup DifferentialEquationSubModule
//! \brief Interface for classes calculating the bounds of a flow.
class BounderInterface {
  public:
    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, Suggestion<StepSizeType> const& hsug) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t)\f$ starting in \f$dom\f$  for time step \f$h\leq h_{sug}\f$.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, Suggestion<StepSizeType> const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! <br>
    //! Arguments: \f$f\f$ is the \a vector_field, \f$D\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& vector_field, BoxDomainType const& state_domain, BoxDomainType const& parameter_domain, Suggestion<StepSizeType> const& suggested_time_step) const = 0;

    //! \brief Compute a pair \f$(h,B)\f$ consisting of a bound \a B for the flow
    //! of \f$\dt{x}=f(x,t,a)\f$ starting in \f$D\f$ with \f$a\in A\f$ for time step \f$h\leq h_{sug}\f$.
    //! Arguments: \f$f\f$ is the \a vector_field, \f$dom\f$ is the \a state_domain and \f$h_{sug}\f$ is the suggested time step.
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& differential_equation, BoxDomainType const& state_domain, StepSizeType const& initial_time, BoxDomainType const& parameter_domain, Suggestion<StepSizeType> const& suggested_time_step) const = 0;

  public:
    /** @brief Perform the integration operation. */
    virtual Void _write(OutputStream& os) const = 0;
    /** @brief Perform the integration operation. */
    virtual BounderInterface* clone() const = 0;

    /** @brief Apply the corresponding integration operator. */
    friend inline OutputStream& operator<<(OutputStream& os, BounderInterface const& bounder) {
        bounder._write(os); return os; }

    virtual ~BounderInterface() = default;
};

/** @brief Type representing bounder base in integration algorithms. */
class BounderBase : public BounderInterface {
  public:
    /** @brief Construct the integration object. */
    BounderBase(LipschitzTolerance lipschitz, MinimumStepSize minimum_step);
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, Suggestion<StepSizeType> const& hsug) const override;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, Suggestion<StepSizeType> const& hsug) const override;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const override = 0;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const override = 0;
  protected:
    const LipschitzTolerance _lipschitz_tolerance;
    const MinimumStepSize _minimum_step_size;
};

//! \ingroup DifferentialEquationSubModule
//! \brief Compute bounds on the flow of a differential equation using a set-based Euler method.
class EulerBounder final : public BounderBase {
  public:
    /** @brief Construct the integration object. */
    EulerBounder(LipschitzTolerance lipschitz = DEFAULT_LIPSCHITZ_TOLERANCE, MinimumStepSize = DEFAULT_MINIMUM_STEP_SIZE);
    using BounderBase::compute;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const override;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const override;
    virtual Void _write(OutputStream& os) const override { os << "EulerBounder"; }
    virtual EulerBounder* clone() const override { return new EulerBounder(*this); }
  private:
    // Compute the bounds on the reach set of dx/dt = f(x,t,a) starting at time t, for a suggested step size of h.
    /** @brief Perform the integration operation. */
    Pair<StepSizeType,UpperBoxType> _compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const;
    // Compute a bounding box for the flow of f starting in D at T.lower() to time T.upper() with parameters A, assuming flow remains in B,
    // but widening D by INITIAL_WIDENING and the flow vectors by VECTOR_WIDENING.
    /** @brief Perform the integration operation. */
    UpperBoxType _formula(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B, PositiveFloatDP INITIAL_WIDENING, PositiveFloatDP VECTOR_WIDENING) const;
    /** @brief Perform the integration operation. */
    UpperBoxType _refinement(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, IntervalDomainType const& T, BoxDomainType const& A, UpperBoxType const& B) const;
};

/** @brief Type representing bounder handle in integration algorithms. */
class BounderHandle {
  private:
    SharedPointer<BounderInterface> _impl;
  public:
    BounderHandle(BounderInterface const& bounder) : _impl(bounder.clone()) { }
    BounderHandle(BounderHandle const& other) : _impl(other._impl) { }
    BounderHandle& operator=(BounderHandle const& other) { _impl = other._impl; return *this; }

    operator BounderInterface const& () const { return *_impl; }

    Void _write(OutputStream& os) const { _impl->_write(os); }
    BounderHandle* clone() const { return new BounderHandle(*this); }

    /** @brief Compute the requested flow, enclosure, or integration step. */
    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, Suggestion<StepSizeType> const& hsug) const {
        return _impl->compute(f,D,hsug);
    }

    /** @brief Compute the requested flow, enclosure, or integration step. */
    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const {
        return _impl->compute(f,D,A,hsug);
    }

    /** @brief Compute the requested flow, enclosure, or integration step. */
    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, Suggestion<StepSizeType> const& hsug) const {
        return _impl->compute(f,D,t,hsug);
    }

    /** @brief Compute the requested flow, enclosure, or integration step. */
    Pair<StepSizeType,UpperBoxType> compute(ValidatedVectorMultivariateFunction const& f, BoxDomainType const& D, StepSizeType const& t, BoxDomainType const& A, Suggestion<StepSizeType> const& hsug) const {
        return _impl->compute(f,D,t,A,hsug);
    }

};



} // namespace Ariadne

#endif /* ARIADNE_BOUNDER_HPP */
