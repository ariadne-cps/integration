/***************************************************************************
 *            test_integrator.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <sstream>
#include <memory>
#include <array>
#include "integration/inclusion_integrator.hpp"
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.tpl.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
List<InputApproximation> approximations() {
    return {InputApproximation(ZeroApproximation()),InputApproximation(ConstantApproximation()),
            InputApproximation(AffineApproximation()),InputApproximation(SinusoidalApproximation()),
            InputApproximation(PiecewiseApproximation())};
}

Void test_approximation_metadata() {
    auto x=EffectiveFormula::coordinate(0u);
    auto u=EffectiveFormula::coordinate(1u);
    EffectiveVectorMultivariateFunction field=EffectiveVectorMultivariateFormulaFunction(2u,Vector<EffectiveFormula>{x+u});
    BoxDomainType inputs={{-0.125_x,0.125_x}};
    TaylorSeriesIntegrator base(step_maximum_error=1e-6,order=8);
    InclusionIntegratorFactory factory(base);
    std::unique_ptr<InclusionIntegratorFactory> cloned_factory(factory.clone());
    List<InclusionIntegrator> integrators;
    Nat index=0u;
    for(auto const& approximation:approximations()) {
        std::unique_ptr<InputApproximationInterface> owned_approximation(approximation.reference().clone());
        ARIADNE_TEST_EQUALS(owned_approximation->index(),index)
        ARIADNE_TEST_EQUALS(owned_approximation->is_reach_rigorous(),index<2u)
        auto integrator=cloned_factory->create(field,inputs,approximation);
        ARIADNE_TEST_EQUALS(integrator.index(),index)
        ARIADNE_TEST_EQUALS(integrator.num_params_per_input(),index==0u ? 0u : (index==1u ? 1u : 2u))
        ARIADNE_TEST_EQUALS(integrator.is_reach_rigorous(),index<2u)
        ARIADNE_TEST_EQUALS(approximation.is_reach_rigorous(),index<2u)
        std::ostringstream expected,actual;
        expected << approximation;
        actual << integrator;
        ARIADNE_TEST_EQUALS(actual.str(),expected.str())
        std::unique_ptr<InclusionIntegrator> cloned(integrator.clone());
        ARIADNE_TEST_ASSERT(*cloned==integrator)
        InclusionIntegratorInterface const& interface=*cloned;
        ARIADNE_TEST_EQUALS(interface.index(),index)
        integrators.append(integrator);
        ++index;
    }
    for(SizeType i=0u;i<integrators.size();++i) {
        for(SizeType j=0u;j<integrators.size();++j) {
            ARIADNE_TEST_EQUALS(integrators[i]==integrators[j],i==j)
            ARIADNE_TEST_EQUALS(integrators[i]<integrators[j],i<j)
        }
    }
    for(auto relation:{InputsRelationKind::AFFINE,InputsRelationKind::SINGULAR,
                       InputsRelationKind::DUAL,InputsRelationKind::ADDITIVE}) {
        std::ostringstream text;
        text << relation;
        ARIADNE_TEST_ASSERT(not text.str().empty())
    }
    ARIADNE_TEST_THROWS((std::ostringstream{} << static_cast<InputsRelationKind>(99)),std::runtime_error)
}

template<class A> Void check_error_processor(EffectiveVectorMultivariateFunction const& f, BoxDomainType const& inputs) {
    ApproximationErrorProcessorFactory<A> factory;
    auto processor=factory.create(f,inputs);
    UpperBoxType bound={{0.25_x,0.5_x},{0.25_x,0.5_x}};
    auto errors=processor->process(PositiveFloatDP(0.0078125_x,dp),bound);
    auto smaller=processor->process(PositiveFloatDP(0.00390625_x,dp),bound);
    ARIADNE_TEST_EQUALS(errors.size(),2u)
    for(SizeType j=0u;j<errors.size();++j) {
        ARIADNE_TEST_ASSERT(std::isfinite(errors[j].raw().get_d()))
        ARIADNE_TEST_ASSERT(errors[j].raw()>=FloatDP(0,dp))
        ARIADNE_TEST_ASSERT(smaller[j].raw()<=errors[j].raw())
    }
}

Void test_error_processor_relations() {
    // Each relation selects a different mathematical error bound. Nonlinear
    // state dependence exercises first and second state derivatives.
    for(SizeType m:{1u,2u,3u}) {
        ARIADNE_PRINT_TEST_CASE_TITLE("State-dependent noise with " << m << " inputs")
        auto x=EffectiveFormula::coordinate(0u);
        auto y=EffectiveFormula::coordinate(1u);
        auto f0=x+x*x+x*y;
        auto f1=y+x*y+y*y;
        BoxDomainType inputs(m,ExactIntervalType(-0.125_x,0.125_x));
        for(SizeType i=0u;i<m;++i) {
            auto u=EffectiveFormula::coordinate(2u+i);
            f0=f0+(1+x+y*y)*u;
            f1=f1+(1+y+x*x)*u;
        }
        EffectiveVectorMultivariateFunction f=EffectiveVectorMultivariateFormulaFunction(2u+m,Vector<EffectiveFormula>{f0,f1});
        check_error_processor<ZeroApproximation>(f,inputs);
        check_error_processor<ConstantApproximation>(f,inputs);
        check_error_processor<AffineApproximation>(f,inputs);
        check_error_processor<SinusoidalApproximation>(f,inputs);
        check_error_processor<PiecewiseApproximation>(f,inputs);
    }
    auto x=EffectiveFormula::coordinate(0u);
    auto y=EffectiveFormula::coordinate(1u);
    auto u=EffectiveFormula::coordinate(2u);
    auto v=EffectiveFormula::coordinate(3u);
    EffectiveVectorMultivariateFunction additive=EffectiveVectorMultivariateFormulaFunction(4u,Vector<EffectiveFormula>{x+x*x+x*y+u,y+x*y+y*y+v});
    BoxDomainType inputs={{-0.125_x,0.125_x},{-0.125_x,0.125_x}};
    check_error_processor<ZeroApproximation>(additive,inputs);
    check_error_processor<ConstantApproximation>(additive,inputs);
    check_error_processor<AffineApproximation>(additive,inputs);
    check_error_processor<SinusoidalApproximation>(additive,inputs);
    check_error_processor<PiecewiseApproximation>(additive,inputs);
}

Void test_nominal_solution_in_reach_and_evolve() {
    auto x=EffectiveFormula::coordinate(0u);
    auto u=EffectiveFormula::coordinate(1u);
    EffectiveVectorMultivariateFunction field=EffectiveVectorMultivariateFormulaFunction(2u,Vector<EffectiveFormula>{x+u});
    BoxDomainType inputs={{-0.125_x,0.125_x}};
    BoxDomainType domain={{0.0_x,0.125_x}};
    UpperBoxType bound={{-1.0_x,1.0_x}};
    // Verify enclosure semantics independently of an optional accuracy budget.
    TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp,1e-10),lipschitz_tolerance=0.5_x,order=8);
    InclusionIntegratorFactory factory(base);
    auto initial=base.function_factory().create_identity(domain);
    StepSizeType h(0.03125_x);
    auto expected=FloatDPBounds(0.0625_x,dp)*exp(FloatDPBounds(h,dp));
    UpperIntervalType expected_interval(expected.lower(),expected.upper());
    for(ExactDouble amplitude:{0.125_x,0.0_x}) {
      BoxDomainType current_inputs={{-amplitude,amplitude}};
      for(auto const& approximation:approximations()) {
        ARIADNE_PRINT_TEST_CASE_TITLE("Nominal trajectory under " << approximation << ", amplitude " << amplitude)
        auto integrator=factory.create(field,current_inputs,approximation);
        auto reaches=integrator.reach(domain,initial,bound,TimeStepType(0),h);
        ARIADNE_TEST_ASSERT(not reaches.empty())
        for(auto const& reach:reaches) {
            ARIADNE_TEST_EQUALS(reach.result_size(),1u)
            ARIADNE_TEST_ASSERT(definitely(is_bounded(reach.range())))
            auto evolved=integrator.evolve(reach,h);
            ARIADNE_TEST_EQUALS(evolved.result_size(),1u)
            ARIADNE_TEST_ASSERT(refines(expected_interval,evolved.range()[0u]))
        }
      }
    }
}

class UnsupportedApproximation : public InputApproximationInterface {
  public:
    bool is_reach_rigorous() const override { return false; }
    Void write(OutputStream& os) const override { os << "UNSUPPORTED"; }
    UnsupportedApproximation* clone() const override { return new UnsupportedApproximation(*this); }
    Nat index() const override { return 99u; }
};

Void test_unsupported_approximation() {
    auto x=EffectiveFormula::coordinate(0u),u=EffectiveFormula::coordinate(1u);
    EffectiveVectorMultivariateFunction field=EffectiveVectorMultivariateFormulaFunction(2u,Vector<EffectiveFormula>{x+u});
    BoxDomainType inputs={{-0.125_x,0.125_x}};
    TaylorSeriesIntegrator base(step_maximum_error=1e-6,order=8);
    InclusionIntegratorFactory factory(base);
    ARIADNE_TEST_THROWS(factory.create(field,inputs,InputApproximation(UnsupportedApproximation())),std::runtime_error)
}

ErrorConstants make_constants(std::array<Vector<ErrorType>,12u> const& v) {
    ErrorType z(0u,dp);
    return ErrorConstants(z,v[0u],z,z,z,v[1u],v[2u],v[3u],
                          z,v[4u],z,z,z,v[5u],v[6u],v[7u],
                          z,v[8u],z,z,z,v[9u],v[10u],v[11u],
                          FloatDPUpperBound(1,dp),ErrorType(2u,dp),ErrorType(3u,dp));
}

Void test_error_constant_dimensions_and_values() {
    std::array<Vector<ErrorType>,12u> vectors;
    for(auto& v:vectors) { v=Vector<ErrorType>(2u,ErrorType(0u,dp)); }
    auto constants=make_constants(vectors);
    ARIADNE_TEST_EQUALS(constants.dimension(),2u)
    auto values=constants.values();
    ARIADNE_TEST_EQUALS(std::get<1u>(values).size(),2u)
    ARIADNE_TEST_EQUALS(std::get<25u>(values).raw(),FloatDP(2,dp))
    ARIADNE_TEST_EQUALS(std::get<26u>(values).raw(),FloatDP(3,dp))
    std::ostringstream text;
    text << constants;
    ARIADNE_TEST_ASSERT(text.str().find("Lambda=")!=std::string::npos)
    // Kj establishes the dimension; each of the remaining eleven vectors
    // must agree, and a caller can supply each mismatch independently.
    for(SizeType j=1u;j<vectors.size();++j) {
        auto invalid=vectors;
        invalid[j]=Vector<ErrorType>(1u,ErrorType(0u,dp));
        ARIADNE_TEST_THROWS(make_constants(invalid),std::runtime_error)
    }
}

Void test_error_dimension_contracts() {
    auto x=EffectiveFormula::coordinate(0u),u=EffectiveFormula::coordinate(1u);
    EffectiveVectorMultivariateFunction field=EffectiveVectorMultivariateFormulaFunction(2u,Vector<EffectiveFormula>{x+u});
    BoxDomainType wrong_inputs(0u);
    ARIADNE_TEST_THROWS(ApproximationErrorProcessorFactory<ZeroApproximation>().create(field,wrong_inputs),std::runtime_error)
    ARIADNE_TEST_THROWS(ApproximationErrorProcessorFactory<ConstantApproximation>().create(field,wrong_inputs),std::runtime_error)
    ARIADNE_TEST_THROWS(ApproximationErrorProcessorFactory<AffineApproximation>().create(field,wrong_inputs),std::runtime_error)
    ARIADNE_TEST_THROWS(ApproximationErrorProcessorFactory<SinusoidalApproximation>().create(field,wrong_inputs),std::runtime_error)
    ARIADNE_TEST_THROWS(ApproximationErrorProcessorFactory<PiecewiseApproximation>().create(field,wrong_inputs),std::runtime_error)
    TaylorSeriesIntegrator base(step_maximum_error=1e-6,order=8);
    BoxDomainType domain={{0.0_x,0.125_x}};
    auto patch=base.function_factory().create_identity(domain);
    Vector<ErrorType> wrong_errors(2u,ErrorType(0u,dp));
    ARIADNE_TEST_THROWS(add_errors(patch,wrong_errors),std::runtime_error)
    ARIADNE_TEST_EQUALS(patch.result_size(),1u)
}

Void test_contracting_drift_zero_approximation() {
    auto x=EffectiveFormula::coordinate(0u),u=EffectiveFormula::coordinate(1u);
    EffectiveVectorMultivariateFunction field=EffectiveVectorMultivariateFormulaFunction(2u,Vector<EffectiveFormula>{-x+u});
    BoxDomainType inputs={{-0.125_x,0.125_x}};
    UpperBoxType bound={{-1.0_x,1.0_x}};
    auto processor=ApproximationErrorProcessorFactory<ZeroApproximation>().create(field,inputs);
    // Additive inputs are also affine. The general processor must preserve
    // the result when this special case is selected explicitly by a caller.
    ApproximationErrorProcessor<ZeroApproximation,AffineInputs> general(field,inputs);
    PositiveFloatDP h(0.03125_x,dp);
    auto errors=processor->process(h,bound);
    auto general_errors=general.process(h,bound);
    ARIADNE_TEST_EQUALS(general_errors[0u].raw(),errors[0u].raw())
    auto exact_error=FloatDPBounds(0.125_x,dp)*(1-exp(-FloatDPBounds(h)));
    ARIADNE_TEST_EQUALS(errors.size(),1u)
    ARIADNE_TEST_ASSERT(errors[0u].raw()>=exact_error.upper().raw())
}

Void test_numeric_helpers() {
    ARIADNE_TEST_EQUALS(lower_bound(Real(0.125_x)),TimeStepType(0.125_x))
    Vector<ErrorType> errors({ErrorType(1u,dp),ErrorType(2u,dp)});
    auto const& raw=cast_exact(errors);
    ARIADNE_TEST_EQUALS(raw.size(),2u)
    ARIADNE_TEST_EQUALS(raw[0u],FloatDP(1,dp))
    ARIADNE_TEST_EQUALS(raw[1u],FloatDP(2,dp))
    errors[0u]=ErrorType(4u,dp);
    ARIADNE_TEST_EQUALS(raw[0u],FloatDP(1,dp))
    Vector<ErrorType> empty_errors(0u,dp);
    ARIADNE_TEST_EQUALS(cast_exact(empty_errors).size(),0u)
    auto const& from_temporary=cast_exact(Vector<ErrorType>({ErrorType(3u,dp)}));
    ARIADNE_TEST_EQUALS(from_temporary[0u],FloatDP(3,dp))
    for(Int value:{-1,1}) {
        FloatDPBounds point(value,dp);
        auto expected=(exp(point)-1)/point;
        auto computed=dexp(point);
        ARIADNE_TEST_ASSERT(refines(expected,computed))
    }
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(1,dp),dexp(FloatDPBounds(0,dp))))
    FloatDPBounds crossing(FloatDP(-1,dp),FloatDP(1,dp));
    auto range=dexp(crossing);
    ARIADNE_TEST_ASSERT(refines(dexp(FloatDPBounds(-1,dp)),range))
    ARIADNE_TEST_ASSERT(refines(dexp(FloatDPBounds(1,dp)),range))
}
}

int main() {
    ARIADNE_TEST_CALL(test_approximation_metadata())
    ARIADNE_TEST_CALL(test_numeric_helpers())
    ARIADNE_TEST_CALL(test_error_constant_dimensions_and_values())
    ARIADNE_TEST_CALL(test_error_dimension_contracts())
    ARIADNE_TEST_CALL(test_contracting_drift_zero_approximation())
    ARIADNE_TEST_CALL(test_unsupported_approximation())
    for(unsigned verbosity:{0u,3u}) {
        Logger::instance().configuration().set_verbosity(verbosity);
        ARIADNE_TEST_CALL(test_error_processor_relations())
        ARIADNE_TEST_CALL(test_nominal_solution_in_reach_and_evolve())
    }
    return ARIADNE_TEST_FAILURES;
}
