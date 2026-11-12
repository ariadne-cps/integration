/***************************************************************************
 *            integration/inclusion_integrator.hpp
 *
 *  Copyright  2008-18  Luca Geretti, Pieter Collins, Sanja Zivanovic
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

/** @file integration/inclusion_integrator.hpp
 *  @brief Input approximations and validated integrators for differential inclusions.
 */

#ifndef ARIADNE_INCLUSION_INTEGRATOR_HPP
#define ARIADNE_INCLUSION_INTEGRATOR_HPP

#include "algebra/algebra.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/vector.hpp"
#include "function/domain.hpp"
#include "function/formula.hpp"
#include "function/function_model.hpp"
#include "function/function_patch.hpp"
#include "function/symbolic_function.hpp"
#include "integration/integrator_interface.hpp"
#include "logging/logging.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

/** @brief Type representing real in integration algorithms. */
class Real;

using TimeStepType = Dyadic;

typedef FloatDPError ErrorType;

/** @brief Perform the integration operation. */
template<class C, class T>
Bool instance_of(T* const obj) {
    return (dynamic_cast<const C*>(obj) != nullptr);
}

/** @brief Perform the integration operation. */
inline TimeStepType lower_bound(const TimeStepType& t) {
    return TimeStepType(FloatDPLowerBound(t, DoublePrecision()).raw());
}
/** @brief Perform the integration operation. */
inline TimeStepType lower_bound(const Real& t) {
    return TimeStepType(FloatDPLowerBound(t, DoublePrecision()).raw());
}

/** @brief Perform the integration operation. */
inline Vector<FloatDP> cast_exact(const Vector<ErrorType>& v) {
    // Vector specializations need not share a layout. Extract the exact values
    // into independent storage instead of aliasing the error vector.
    Vector<FloatDP> result(v.size(), v.zero_element().raw());
    for (SizeType i = 0u; i < v.size(); ++i) {
        result[i] = v[i].raw();
    }
    return result;
}

/** @brief Perform the integration operation. */
Void add_errors(ValidatedVectorMultivariateFunctionPatch& phi, const Vector<ErrorType>& e);

/** @brief Perform the integration operation. */
EffectiveVectorMultivariateFunction substitute_v_with_w(const EffectiveVectorMultivariateFunction& F, const Vector<EffectiveScalarMultivariateFunction>& w);

template<class F1, class F2, class F3, class... FS>
decltype(auto) product(F1 const& f1, F2 const& f2, F3 const& f3, FS const&... fs) {
    return product(product(f1, f2), f3, fs...);
}

template<class F1, class F2, class F3, class... FS>
decltype(auto) join(F1 const& f1, F2 const& f2, F3 const& f3, FS const&... fs) {
    return join(join(f1, f2), f3, fs...);
}

//! \brief The continuous function dexp(x)=(exp(x)-1)/x, with dexp(0)=1.
//! Note that the function is positive and monotone increasing.
template<class F>
PositiveUpperBound<F> dexp(const UpperBound<F>& x) {
    if (x.raw() == 0) {
        return PositiveUpperBound<F>(1u, x.precision());
    }
    // Evaluate the exact endpoint with interval arithmetic. Negation reverses
    // rounding directions, so reusing an upper exp bound at negative x is unsafe.
    /** @brief Perform the integration operation. */
    const Bounds<F> point(x.raw());
    return cast_positive((exp(point) - 1) / point).upper();
}

/** @brief Perform the integration operation. */
template<class F>
PositiveLowerBound<F> dexp(const LowerBound<F>& x) {
    if (x.raw() == 0) {
        return PositiveLowerBound<F>(1u, x.precision());
    }
    /** @brief Perform the integration operation. */
    const Bounds<F> point(x.raw());
    return cast_positive((exp(point) - 1) / point).lower();
}

/** @brief Perform the integration operation. */
template<class F>
PositiveBounds<F> dexp(const Bounds<F>& x) {
    return PositiveBounds<F>(dexp(x.lower()), dexp(x.upper()));
}

// Compute (x*e^x/2-e^x+x/2+1)/x^3, which is a positive monotone function.
/** @brief Perform the integration operation. */
template<class F>
PositiveBounds<F> psi0(const F& x) {
    return PositiveBounds<F>((x * exp(x) / 2u - exp(x) + x / 2u + 1u) / pow(x, 3u));
}
// Compute (x*e^x-e^x-x^2/2+1)/x^3, which is a positive monotone function.
/** @brief Perform the integration operation. */
template<class F>
PositiveBounds<F> psi1(const F& x) {
    return PositiveBounds<F>(((x - 1) * exp(x) + (1u - sqr(x) / 2u)) / pow(x, 3u));
}

/** @brief Perform the integration operation. */
template<class F>
PositiveUpperBound<F> psi0(const UpperBound<F>& x) {
    return psi0(cast_exact(x));
}

/** @brief Perform the integration operation. */
template<class F>
PositiveUpperBound<F> psi1(const UpperBound<F>& x) {
    return psi1(cast_exact(x));
}

/** @brief Perform the integration operation. */
template<class F>
PositiveUpperBound<F> psi0(const PositiveUpperBound<F>& xu) {
    PositiveLowerBound<F> xl = cast_exact(xu);
    return cast_positive(xu * exp(xu) / 2u + xu / 2u + 1u - exp(xl)) / pow(xl, 3u);
}

/** @brief Perform the integration operation. */
template<class F>
PositiveUpperBound<F> psi1(const PositiveUpperBound<F>& xu) {
    PositiveLowerBound<F> xl = cast_exact(xu);
    return cast_positive(xu * exp(xu) + 1u - (exp(xl) + sqr(xl) / 2u)) / pow(xl, 3u);
}

/** @brief Constants controlling integration error estimates. */
struct ErrorConstants {
    typedef decltype(declval<ErrorType>().precision()) PrecisionType;

    ErrorType K;              // K
    Vector<ErrorType> Kj;     // K[j]
    ErrorType pK;             // K'
    ErrorType pKv;            // K'v
    ErrorType pKw;            // K'w
    Vector<ErrorType> pKj;    // K'[j]
    Vector<ErrorType> pKjv;   // K'v[j]
    Vector<ErrorType> pKjw;   // K'w[j]
    ErrorType L;              // L
    Vector<ErrorType> Lj;     // L[j]
    ErrorType pL;             // L'
    ErrorType pLv;            // L'v
    ErrorType pLw;            // L'w
    Vector<ErrorType> pLj;    // L'[j]
    Vector<ErrorType> pLjv;   // L'v[j]
    Vector<ErrorType> pLjw;   // L'w[j]
    ErrorType H;              // H
    Vector<ErrorType> Hj;     // H[j]
    ErrorType pH;             // H'
    ErrorType pHv;            // H'v
    ErrorType pHw;            // H'w
    Vector<ErrorType> pHj;    // H'[j]
    Vector<ErrorType> pHjv;   // H'v[j]
    Vector<ErrorType> pHjw;   // H'w[j]
    FloatDPUpperBound Lambda; // Lambda
    ErrorType expLambda;      // e^(Lambda*h - 1) / (Lambda*h)
    ErrorType expL;           // e^(L*h)

    ErrorConstants(const ErrorType&, const Vector<ErrorType>&, const ErrorType&, const ErrorType&, const ErrorType&, const Vector<ErrorType>&, const Vector<ErrorType>&, const Vector<ErrorType>&,
                   const ErrorType&, const Vector<ErrorType>&, const ErrorType&, const ErrorType&, const ErrorType&, const Vector<ErrorType>&, const Vector<ErrorType>&, const Vector<ErrorType>&,
                   const ErrorType&, const Vector<ErrorType>&, const ErrorType&, const ErrorType&, const ErrorType&, const Vector<ErrorType>&, const Vector<ErrorType>&, const Vector<ErrorType>&, const FloatDPUpperBound&, const ErrorType&, const ErrorType&);
    Tuple<ErrorType, Vector<ErrorType>, ErrorType, ErrorType, ErrorType, Vector<ErrorType>, Vector<ErrorType>, Vector<ErrorType>,
          ErrorType, Vector<ErrorType>, ErrorType, ErrorType, ErrorType, Vector<ErrorType>, Vector<ErrorType>, Vector<ErrorType>,
          ErrorType, Vector<ErrorType>, ErrorType, ErrorType, ErrorType, Vector<ErrorType>, Vector<ErrorType>, Vector<ErrorType>,
          /** @brief Return the requested integration parameter or result. */
          FloatDPUpperBound, ErrorType, ErrorType>
        values() const;

    SizeType dimension() const { return _dimension; }
    PrecisionType precision() const { return PrecisionType(); }

  private:
    SizeType _dimension;
};

inline std::ostream& operator<<(std::ostream& os, const ErrorConstants& n) {
    os << "K=" << n.K << ", Kj=" << n.Kj << ", K'=" << n.pK << ", K'v=" << n.pKv << ", Kw'=" << n.pKw << ", Kj'=" << n.Kj << ", K'jv=" << n.pKjv << ", K'jw=" << n.pKjw << ", L=" << n.L << ", Lj=" << n.Lj << ", L'=" << n.pL << ", L'v=" << n.pLv << ", L'w=" << n.pLw << ", Lj'=" << n.pLj << ", Lj'v=" << n.pLjv << ", Lj'w=" << n.pLjw << ", H=" << n.H << ", Hj=" << n.Hj << ", H'=" << n.pH << ", H'v=" << n.pHv << ", H'w=" << n.pHw << ", Hj'=" << n.pHj << ", Hj'v=" << n.pHjv << ", Hj'w=" << n.pHjw << ", Lambda=" << n.Lambda << ", expLambda=" << n.expLambda << ", expL=" << n.expL;
    return os;
}

/** @brief Compute the requested flow, enclosure, or integration step. */
ErrorConstants compute_constants(const EffectiveVectorMultivariateFunction&, const Vector<EffectiveVectorMultivariateFunction>&, const BoxDomainType&, const PositiveFloatDP&, const UpperBoxType&);

/** @brief Abstract interface for the temporal approximation of an input. */
struct InputApproximationInterface {
    //! \brief Whether the reached set is rigorous or another approximation must be used
    virtual bool is_reach_rigorous() const = 0;
    /** @brief Write an integration object or its results to an output stream. */
    virtual Void write(OutputStream& os) const = 0;
    /** @brief Perform the integration operation. */
    virtual InputApproximationInterface* clone() const = 0;
    virtual ~InputApproximationInterface() = default;
    /** @brief Perform the integration operation. */
    virtual Nat index() const = 0;
};

inline std::ostream& operator<<(std::ostream& os, const InputApproximationInterface& a) {
    a.write(os);
    return os;
}

/** @brief Input approximation identically equal to zero. */
struct ZeroApproximation : public InputApproximationInterface {
    virtual bool is_reach_rigorous() const override { return true; }
    virtual Void write(OutputStream& os) const override { os << "ZERO"; }
    virtual Nat index() const override { return 0; }
    virtual ZeroApproximation* clone() const override { return new ZeroApproximation(*this); }
};
/** @brief Constant-valued input approximation. */
struct ConstantApproximation : public InputApproximationInterface {
    virtual bool is_reach_rigorous() const override { return true; }
    virtual Void write(OutputStream& os) const override { os << "CONSTANT"; }
    virtual Nat index() const override { return 1; }
    virtual ConstantApproximation* clone() const override { return new ConstantApproximation(*this); }
};
/** @brief Affine-in-time input approximation. */
struct AffineApproximation : public InputApproximationInterface {
    virtual bool is_reach_rigorous() const override { return false; }
    virtual Void write(OutputStream& os) const override { os << "AFFINE"; }
    virtual Nat index() const override { return 2; }
    virtual AffineApproximation* clone() const override { return new AffineApproximation(*this); }
};
/** @brief Sinusoidal input approximation. */
struct SinusoidalApproximation : public InputApproximationInterface {
    virtual bool is_reach_rigorous() const override { return false; }
    virtual Void write(OutputStream& os) const override { os << "SINUSOIDAL"; }
    virtual Nat index() const override { return 3; }
    virtual SinusoidalApproximation* clone() const override { return new SinusoidalApproximation(*this); }
};
/** @brief Piecewise-defined input approximation. */
struct PiecewiseApproximation : public InputApproximationInterface {
    virtual bool is_reach_rigorous() const override { return false; }
    virtual Void write(OutputStream& os) const override { os << "PIECEWISE"; }
    virtual Nat index() const override { return 4; }
    virtual PiecewiseApproximation* clone() const override { return new PiecewiseApproximation(*this); }
};

/** @brief Value-like handle for a chosen input approximation. */
class InputApproximation : public Handle<const InputApproximationInterface> {
  public:
    using Handle<const InputApproximationInterface>::Handle;
    bool is_reach_rigorous() const { return this->_ptr->is_reach_rigorous(); }
    template<class A>
    Bool handles(const A&) const { return instance_of<A>(&*this->_ptr); }
};

/** @brief Perform the integration operation. */
template<class A>
constexpr Nat const_num_params_per_input();
template<>
constexpr Nat const_num_params_per_input<ZeroApproximation>() { return 0u; }
template<>
constexpr Nat const_num_params_per_input<ConstantApproximation>() { return 1u; }
template<>
constexpr Nat const_num_params_per_input<AffineApproximation>() { return 2u; }
template<>
constexpr Nat const_num_params_per_input<SinusoidalApproximation>() { return 2u; }
template<>
constexpr Nat const_num_params_per_input<PiecewiseApproximation>() { return 2u; }

/** @brief Relationship between the differential inclusion inputs. */
enum class InputsRelationKind : std::uint8_t { AFFINE,
                                               SINGULAR,
                                               DUAL,
                                               ADDITIVE };

template<InputsRelationKind R>
struct InputsRelationKindTrait {
    static InputsRelationKind kind() { return R; }
};

/** @brief Perform the integration operation. */
ErrorType vstar(const BoxDomainType& inputs);

/** @brief Perform the integration operation. */
template<class A>
ErrorType wstar_multiplier();

/** @brief Perform the integration operation. */
template<class A>
ErrorType wstar(const BoxDomainType& inputs) {
    typedef typename ErrorType::PrecisionType PR;
    PR pr;
    ErrorType result(0u, pr);
    for (auto i: range(0, inputs.size()))
        result = max(result, wstar_multiplier<A>() * ErrorType(abs(inputs[i]).upper_bound()));
    return result;
}

/** @brief Tag identifying affine input relationships. */
class AffineInputs : public InputsRelationKindTrait<InputsRelationKind::AFFINE> {};
/** @brief Tag identifying a singular input relationship. */
class SingularInput : public InputsRelationKindTrait<InputsRelationKind::SINGULAR> {};
/** @brief Tag identifying dual input relationships. */
class DualInputs : public InputsRelationKindTrait<InputsRelationKind::DUAL> {};
/** @brief Tag identifying additive input relationships. */
class AdditiveInputs : public InputsRelationKindTrait<InputsRelationKind::ADDITIVE> {};

/** @brief Apply the corresponding integration operator. */
inline std::ostream& operator<<(std::ostream& os, const InputsRelationKind& kind) {
    switch (kind) {
        case InputsRelationKind::AFFINE:
            os << "AFFINE";
            break;
        case InputsRelationKind::SINGULAR:
            os << "SINGULAR";
            break;
        case InputsRelationKind::DUAL:
            os << "DUAL";
            break;
        case InputsRelationKind::ADDITIVE:
            os << "ADDITIVE";
            break;
        /** @brief Construct the integration object. */
        default:
            ARIADNE_FAIL_MSG("Unhandled InputsRoles for output streaming");
    }
    return os;
}

/** @brief Perform the integration operation. */
template<class A>
Vector<EffectiveScalarMultivariateFunction> build_w_functions(const Interval<TimeStepType>& domt, const BoxDomainType& doma, SizeType n, SizeType m);

/** @brief Handle for differential-inclusion integration. */
class InclusionIntegrator;
/** @brief Interface for differential-inclusion evolution. */
class InclusionIntegratorInterface;

/** @brief Factory for differential-inclusion integration strategies. */
class InclusionIntegratorFactory {
    SharedPointer<IntegratorInterface> _integrator;

  public:
    InclusionIntegratorFactory(const IntegratorInterface& integrator) : _integrator(integrator.clone()) {}
    /** @brief Perform the integration operation. */
    InclusionIntegrator create(const EffectiveVectorMultivariateFunction& f, const BoxDomainType& inputs, const InputApproximation& approximation) const;
    InclusionIntegratorFactory* clone() const { return new InclusionIntegratorFactory(*this); }
};

template<class A>
class ApproximationErrorProcessorInterface {
  public:
    /** @brief Perform the integration operation. */
    virtual Vector<ErrorType> process(const PositiveFloatDP& h, const UpperBoxType& B) const = 0;
};

template<class A, class R>
/** @brief Processes approximation errors for differential inclusions. */
class ApproximationErrorProcessor : public ApproximationErrorProcessorInterface<A> {
  public:
    ApproximationErrorProcessor(const EffectiveVectorMultivariateFunction& f, const BoxDomainType& inputs)
        : _f(f), _inputs(inputs) {}
    /** @brief Perform the integration operation. */
    virtual Vector<ErrorType> process(const PositiveFloatDP& h, const UpperBoxType& B) const override;

  private:
    const EffectiveVectorMultivariateFunction& _f;
    const BoxDomainType& _inputs;

  private:
    /** @brief Perform the integration operation. */
    Vector<ErrorType> process(const ErrorConstants& n, const PositiveFloatDP& h) const;

  public:
    virtual ~ApproximationErrorProcessor() = default;
};

// These supported implementations are provided by the shared library.
extern template class ApproximationErrorProcessor<ZeroApproximation, AffineInputs>;
extern template class ApproximationErrorProcessor<ZeroApproximation, SingularInput>;
extern template class ApproximationErrorProcessor<ZeroApproximation, DualInputs>;
extern template class ApproximationErrorProcessor<ZeroApproximation, AdditiveInputs>;
extern template class ApproximationErrorProcessor<ConstantApproximation, AffineInputs>;
extern template class ApproximationErrorProcessor<ConstantApproximation, SingularInput>;
extern template class ApproximationErrorProcessor<ConstantApproximation, DualInputs>;
extern template class ApproximationErrorProcessor<ConstantApproximation, AdditiveInputs>;
extern template class ApproximationErrorProcessor<AffineApproximation, AffineInputs>;
extern template class ApproximationErrorProcessor<AffineApproximation, SingularInput>;
extern template class ApproximationErrorProcessor<AffineApproximation, DualInputs>;
extern template class ApproximationErrorProcessor<AffineApproximation, AdditiveInputs>;
extern template class ApproximationErrorProcessor<SinusoidalApproximation, AffineInputs>;
extern template class ApproximationErrorProcessor<SinusoidalApproximation, SingularInput>;
extern template class ApproximationErrorProcessor<SinusoidalApproximation, DualInputs>;
extern template class ApproximationErrorProcessor<SinusoidalApproximation, AdditiveInputs>;
extern template class ApproximationErrorProcessor<PiecewiseApproximation, AffineInputs>;
extern template class ApproximationErrorProcessor<PiecewiseApproximation, SingularInput>;
extern template class ApproximationErrorProcessor<PiecewiseApproximation, DualInputs>;
extern template class ApproximationErrorProcessor<PiecewiseApproximation, AdditiveInputs>;

template<class A>
/** @brief Factory for approximation-error processors. */
class ApproximationErrorProcessorFactory {
    typedef ApproximationErrorProcessorInterface<A> Processor;

  public:
    /** @brief Perform the integration operation. */
    SharedPointer<Processor> create(const EffectiveVectorMultivariateFunction& f, const BoxDomainType& inputs) const {
        SizeType n = f.result_size();
        SizeType m = inputs.size();
        /** @brief Construct the integration object. */
        ARIADNE_PRECONDITION_MSG(f.argument_size() == n + m, "ApproximationErrorProcessorFactory was given an incompatible f argument space in respect to the inputs box");
        Set<Nat> input_idx;
        for (SizeType i: range(n, n + m)) {
            input_idx.insert(static_cast<Nat>(i));
        }
        if (is_additive_in(f, input_idx))
            return SharedPointer<Processor>(new ApproximationErrorProcessor<A, AdditiveInputs>(f, inputs));
        else if (m == 1)
            return SharedPointer<Processor>(new ApproximationErrorProcessor<A, SingularInput>(f, inputs));
        else if (m == 2)
            return SharedPointer<Processor>(new ApproximationErrorProcessor<A, DualInputs>(f, inputs));
        else
            return SharedPointer<Processor>(new ApproximationErrorProcessor<A, AffineInputs>(f, inputs));
    }
};

/** @brief Interface for differential-inclusion evolution. */
class InclusionIntegratorInterface {
  public:
    /** @brief Write an integration object or its results to an output stream. */
    virtual Void write(OutputStream& os) const = 0;
    virtual Bool operator==(const InclusionIntegratorInterface& rhs) const = 0;
    /** @brief Apply the corresponding integration operator. */
    virtual Bool operator<(const InclusionIntegratorInterface& rhs) const = 0;
    /** @brief Perform the integration operation. */
    virtual Nat index() const = 0;
    /** @brief Perform the integration operation. */
    virtual Nat num_params_per_input() const = 0;
    /** @brief Return the requested integration parameter or result. */
    virtual bool is_reach_rigorous() const = 0;
    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual List<ValidatedVectorMultivariateFunctionPatch> reach(const BoxDomainType& D, const ValidatedVectorMultivariateFunctionPatch& evolve_function, const UpperBoxType& B, const TimeStepType& t, const StepSizeType& h) const = 0;
    /** @brief Perform the integration operation. */
    virtual ValidatedVectorMultivariateFunctionPatch evolve(const ValidatedVectorMultivariateFunctionPatch& reach_function, const TimeStepType& t) const = 0;
    friend std::ostream& operator<<(std::ostream& os, const InclusionIntegratorInterface& approximator) {
        approximator.write(os);
        return os;
    }
    virtual ~InclusionIntegratorInterface() = default;
};

template<class A>
/** @brief Concrete differential-inclusion integrator implementation. */
class InclusionIntegratorImpl : public InclusionIntegratorInterface {
    friend class InclusionIntegratorFactory;

  protected:
    const EffectiveVectorMultivariateFunction& _f;
    const BoxDomainType& _inputs;
    SharedPointer<IntegratorInterface> _integrator;
    SharedPointer<ApproximationErrorProcessorInterface<A>> _processor;
    InclusionIntegratorImpl(const EffectiveVectorMultivariateFunction& f, const BoxDomainType& inputs, const SharedPointer<IntegratorInterface>& integrator) : _f(f), _inputs(inputs), _integrator(integrator), _processor(ApproximationErrorProcessorFactory<A>().create(f, inputs)), _num_params_per_input(const_num_params_per_input<A>()) {}

  private:
    const Nat _num_params_per_input;

  public:
    virtual Void write(OutputStream& os) const override { os << A(); }

    /** @brief Compute the requested flow, enclosure, or integration step. */
    virtual List<ValidatedVectorMultivariateFunctionPatch> reach(const BoxDomainType& D, const ValidatedVectorMultivariateFunctionPatch& evolve_function, const UpperBoxType& B, const TimeStepType& t, const StepSizeType& h) const override;
    /** @brief Perform the integration operation. */
    virtual ValidatedVectorMultivariateFunctionPatch evolve(const ValidatedVectorMultivariateFunctionPatch& reach_function, const TimeStepType& t) const override;

    virtual Bool operator==(const InclusionIntegratorInterface& rhs) const override;
    /** @brief Apply the corresponding integration operator. */
    virtual Bool operator<(const InclusionIntegratorInterface& rhs) const override;

    virtual Nat index() const override { return A().index(); }

    virtual bool is_reach_rigorous() const override { return A().is_reach_rigorous(); }

    virtual Nat num_params_per_input() const override { return _num_params_per_input; }

    virtual ~InclusionIntegratorImpl() = default;

  private:
    /** @brief Compute the requested flow, enclosure, or integration step. */
    Vector<ErrorType> compute_errors(const StepSizeType& h, const UpperBoxType& B) const { return _processor->process(PositiveFloatDP(h, DoublePrecision()), B); };
    /** @brief Perform the integration operation. */
    BoxDomainType build_parameter_domain(const BoxDomainType& V) const;
    /** @brief Perform the integration operation. */
    ValidatedVectorMultivariateFunctionPatch build_reach_function(const ValidatedVectorMultivariateFunctionPatch& evolve_function, const ValidatedVectorMultivariateFunctionPatch& Phi, const TimeStepType& t, const TimeStepType& new_t) const;
    /** @brief Perform the integration operation. */
    ValidatedVectorMultivariateFunctionPatch build_secondhalf_piecewise_reach_function(const ValidatedVectorMultivariateFunctionPatch& evolve_function, const ValidatedVectorMultivariateFunctionPatch& Phi, const TimeStepType& t, const TimeStepType& new_t) const;
    /** @brief Perform the integration operation. */
    Vector<EffectiveScalarMultivariateFunction> build_secondhalf_piecewise_w_functions(const Interval<TimeStepType>& domt, const BoxDomainType& doma, SizeType n, SizeType m) const;
};

/** @brief Handle for differential-inclusion integration. */
class InclusionIntegrator : public Handle<const InclusionIntegratorInterface> {
  public:
    using Handle<const InclusionIntegratorInterface>::Handle;
    virtual ~InclusionIntegrator() = default;
    template<class A>
    Bool handles(const A&) const { return instance_of<A>(&*this->_ptr); }

    friend Bool operator==(const InclusionIntegrator& lhs, const InclusionIntegrator& rhs) {
        return lhs._ptr->operator==(*rhs._ptr);
    }
    /** @brief Apply the corresponding integration operator. */
    friend Bool operator<(const InclusionIntegrator& lhs, const InclusionIntegrator& rhs) {
        return lhs._ptr->operator<(*rhs._ptr);
    }

    bool is_reach_rigorous() const { return _ptr->is_reach_rigorous(); }
    Nat index() const { return _ptr->index(); }
    Nat num_params_per_input() const { return _ptr->num_params_per_input(); }

    friend std::ostream& operator<<(std::ostream& os, const InclusionIntegrator& approximator) {
        os << *approximator._ptr;
        return os;
    }

    operator const InclusionIntegratorInterface&() const { return *_ptr; }
    InclusionIntegrator* clone() const { return new InclusionIntegrator(*this); }

    List<ValidatedVectorMultivariateFunctionPatch> reach(const BoxDomainType& D, const ValidatedVectorMultivariateFunctionPatch& evolve_function, const UpperBoxType& B, const TimeStepType& t, const StepSizeType& h) const { return _ptr->reach(D, evolve_function, B, t, h); }
    ValidatedVectorMultivariateFunctionPatch evolve(const ValidatedVectorMultivariateFunctionPatch& reach_function, const TimeStepType& t) const { return _ptr->evolve(reach_function, t); }
};

} // namespace Ariadne

#endif // ARIADNE_INCLUSION_INTEGRATOR_HPP
