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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include "integration/integrator.hpp"
#include "algebra/vector.hpp"
#include "algebra/algebra.hpp"
#include "algebra/polynomial.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "logging/logging_command_line_interface.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

struct UnsafeReal : Real { UnsafeReal(double d) : Real(operator""_x(d)) { } };

class TestIntegrator : public Ariadne::Testing::ContextFixture
{
    typedef Vector<ExactIntervalType> ExactIntervalVectorType;
  private:
    std::unique_ptr<IntegratorInterface> integrator_ptr;
  public:
    TestIntegrator(Ariadne::Testing::Context& context, const IntegratorInterface& i)
        : ContextFixture(context), integrator_ptr(i.clone())
    { }

    

    

    Void test_constant_derivative() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={o*2,o*3};
        ARIADNE_TESTING_NOTE(f);
        ExactBoxType d={{0.0_x,1.0_x},{-0.5_x,1.5_x}};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,suggest(h));
        EffectiveVectorMultivariateFunction expected_flow={x0+2*t,y0+3*t};
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-8_pr);
    }

    Void test_quadratic_flow() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={o,x};
        ExactBoxType d={ExactIntervalType(0.0_x,1.0_x),ExactIntervalType(-0.5_x,1.5_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,suggest(h));
        EffectiveVectorMultivariateFunction expected_flow={x0+t,y0+x0*t+t*t/2};
        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-8_pr);
    }

    Void test_linear() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        EffectiveVectorMultivariateFunction f={x,-y};
        ExactBoxType d={ExactIntervalType(-0.25_x,0.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,suggest(h));
        EffectiveVectorMultivariateFunction expected_flow={x0*(1+t+t*t/2+t*t*t/6+t*t*t*t/24),y0*(1-t+t*t/2-t*t*t/6+t*t*t*t/24)};
        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-4_pr);
    };

    Void test_spiral() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(2,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(3,0);
        EffectiveScalarMultivariateFunction y0=EffectiveScalarMultivariateFunction::coordinate(3,1);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(3,2);

        Real half(0.5_x);
        EffectiveVectorMultivariateFunction f={-half*x-y,x-half*y};
        ExactBoxType d={ExactIntervalType(0.75_x,1.25_x),ExactIntervalType(-0.25_x,0.25_x)};
        StepSizeType h=0.25_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,suggest(h));
        EffectiveVectorMultivariateFunction expected_flow={exp(-half*t)*(x0*cos(t)-y0*sin(t)),exp(-half*t)*(x0*sin(t)+y0*cos(t))};
        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_NOTE(flow-expected_flow);
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-3_pr);
    }

    // Equation: dx/dt = x*(1-x)
    // Solution:  x = 1/(1+(1/x0-1)*exp(-t))
    Void test_logistic() {

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(1,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(2,1);

        EffectiveVectorMultivariateFunction f={x*(o-x)};
        ExactBoxType d={ExactIntervalType(0.25_x,0.5_x)};
        StepSizeType h=0.5_x;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,suggest(h));
        ValidatedVectorMultivariateTaylorFunctionModelDP taylor_flow=dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP&>(flow.reference());
        EffectiveVectorMultivariateFunction expected_flow={1/(1+(1/x0-1)*exp(-t))};
        ARIADNE_TESTING_NOTE(*integrator_ptr);
        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(taylor_flow.errors());
        ARIADNE_TESTING_NOTE(taylor_flow-expected_flow);
        ARIADNE_TESTING_EXPECT_LESS(taylor_flow.error(), 0.01_pr);
        ARIADNE_TESTING_EXPECT_LESS(norm(taylor_flow-expected_flow), 0.014_pr);
    };

    Void test_time_variant() {

        // dx/dt=a*x+b*t+c; x=(x0+b/a*t0+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a
        // dx/dt=b*t+c; x=x0+c*(t-t0)+b/2*(t^2-t0^2)
        // dx/dt=c; x=x0+c*(t-t0)

        EffectiveScalarMultivariateFunction o=EffectiveScalarMultivariateFunction::constant(4,1);
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(4,0);
        EffectiveScalarMultivariateFunction y=EffectiveScalarMultivariateFunction::coordinate(4,1);
        EffectiveScalarMultivariateFunction z=EffectiveScalarMultivariateFunction::coordinate(4,2);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(4,3);

        Rational a=2; Rational b=3; Rational c=5;
        EffectiveVectorMultivariateFunction f={a*x+b*t+c,b*t+c,c*o};
        //EffectiveVectorMultivariateFunction f={a*x+b*t,b*t+c};
        ExactBoxType d={ExactIntervalType(-0.5_x,1.5_x),ExactIntervalType(-0.5_x,2.5_x),ExactIntervalType(0.0_x,1.0_x)};
        StepSizeType t0=3.0_x;
        StepSizeType hsug=0.0625_x;
        Pair<StepSizeType,UpperBoxType> step_bounds = EulerBounder().compute(f,d,t0,ExactBoxType(0u),suggest(hsug));
        StepSizeType h = step_bounds.first;
        UpperBoxType B = step_bounds.second;
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,d,Interval<StepSizeType>(t0,t0+h),ExactBoxType(0u),B);
        EffectiveVectorMultivariateFunction expected_flow={(x+b/a*t0+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a,y+c*(t-t0)+b/2*(t*t-t0*t0),z+c*(t-t0)};
        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-4_pr);
    };

    Void test_time_variant_with_parameters() {
        // Equation: dx/dt=a*x+b*t+c;
        // Solution: x=c*exp(a*t)-b/a*t-b/a^2
        //     where c=(x0+b/a*t0+b/a^2)/exp(a*t0)
        //         since dx/dt=a*c*exp(a*t)-b/a
        //             a*x+b*t=a*c*exp(a*t)-b*t-b/a + b*t
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(5,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(5,1);
        EffectiveScalarMultivariateFunction a=EffectiveScalarMultivariateFunction::coordinate(5,2);
        EffectiveScalarMultivariateFunction b=EffectiveScalarMultivariateFunction::coordinate(5,3);
        EffectiveScalarMultivariateFunction c=EffectiveScalarMultivariateFunction::coordinate(5,4);
        EffectiveVectorMultivariateFunction f={a*x+b*t+c};
        ExactBoxType domx={ExactIntervalType(0.5_x,1.0_x)};
        StepSizeType t0=0.0_x;
        StepSizeType hsug=0.0625_x;
        ExactBoxType domp={ExactIntervalType{2.0_x,2.5_x},ExactIntervalType{-0.5_x,1.0_x},ExactIntervalType{0.5_x,1.0_x}};
        Pair<StepSizeType,UpperBoxType> step_bounds = EulerBounder().compute(f,domx,t0,domp,suggest(hsug));
        StepSizeType h = step_bounds.first;
        UpperBoxType B = step_bounds.second;
        Interval<StepSizeType> domt(t0,t0+h);
        ValidatedVectorMultivariateFunctionPatch flow=integrator_ptr->flow_step(f,domx,domt,domp,B);
        EffectiveVectorMultivariateFunction expected_flow={(x+b/a*(t0)+(b/a+c)/a)*exp(a*(t-t0))-b/a*t-(b/a+c)/a};
        //   The function below should be zero for phi=expected_flow[0]
        // EffectiveScalarMultivariateFunction phi_error=derivative(phi,1)-(a*phi+b*t+c);

        ARIADNE_TESTING_NOTE(f);
        ARIADNE_TESTING_NOTE(flow);
        ARIADNE_TESTING_NOTE(expected_flow);
        ARIADNE_TESTING_NOTE(flow-expected_flow);
        ARIADNE_TESTING_NOTE(flow.errors());
        ARIADNE_TESTING_EXPECT_LESS(norm(flow-expected_flow), 1e-4_pr);
    }

    // Equation: dx/dt=x^2;
    // Solution: x=1/(1/x0-t)=x0/(1-x0*t)
    Void test_step_size() {
        EffectiveScalarMultivariateFunction x=EffectiveScalarMultivariateFunction::coordinate(1,0);
        EffectiveScalarMultivariateFunction x0=EffectiveScalarMultivariateFunction::coordinate(2,0);
        EffectiveScalarMultivariateFunction t=EffectiveScalarMultivariateFunction::coordinate(2,1);

        EffectiveVectorMultivariateFunction f={sqr(x)};
        EffectiveVectorMultivariateFunction phi={rec(rec(x0)-t)};
        ExactBoxType domx={ExactIntervalType(0.875_x,1.125_x)};
        StepSizeType h=0.0009765625_dy;

        try {
            ValidatedVectorMultivariateFunctionPatch flow(integrator_ptr->flow_step(f,domx,h)); ARIADNE_TESTING_NOTE(flow);
            ARIADNE_TESTING_EXPECT_EQUAL(flow.domain()[0],domx[0]);
            ARIADNE_TESTING_EXPECT_EQUAL(flow.domain()[1].upper_bound(),h);
            ARIADNE_TESTING_NOTE(flow.domain());
            ARIADNE_TESTING_EXPECT_EQUAL(flow.domain()[1].lower_bound(),0);
        } catch (IncompleteFlowException const& e) {
            do { std::ostringstream message; message << "Integrator "<<*integrator_ptr<<" cannot compute flow of dx/dt=x^2 starting in [1.0:1.5] up to time 1/1024."; this->ctx().warning(message.str()); } while(false);
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") threw IncompleteFlowException."; this->ctx().warning(message.str()); } while(false);
            ARIADNE_TESTING_EXPECT_LESS(e.computed_model().domain()[1].upper_bound(), h);
        } catch (FlowBoundsException const&) {
            do { std::ostringstream message; message << "Integrator "<<*integrator_ptr<<" cannot compute flow of dx/dt=x^2 starting in [1.0:1.5] up to time 1/1024."; this->ctx().warning(message.str()); } while(false);
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") threw FlowBoundsException."; this->ctx().warning(message.str()); } while(false);
        }

        h=1.0_dy;

        try {
            ValidatedVectorMultivariateFunctionPatch flow(integrator_ptr->flow_step(f,domx,h)); ARIADNE_TESTING_NOTE(flow);
            ARIADNE_TESTING_EXPECT_EQUAL(flow.domain()[1].upper_bound(),h);
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") completed, which should be impossible."; this->ctx().warning(message.str()); } while(false);
        } catch (IncompleteFlowException const& e) {
            ARIADNE_TESTING_NOTE("Integrator::flow_step(f,d,"<<h<<") caught IncompleteFlowException as expected.");
            ARIADNE_TESTING_EXPECT_LESS(e.computed_model().domain()[1].upper_bound(), h);
        } catch (FlowBoundsException const&) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") threw FlowBoundsException."; this->ctx().warning(message.str()); } while(false);
        } catch (FlowTimeStepException const&) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") threw FlowTimeStepException."; this->ctx().warning(message.str()); } while(false);
        } catch (DivideByZeroException const&) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,"<<h<<") threw DivideByZeroException."; this->ctx().warning(message.str()); } while(false);
        }

        try {
            StepSizeType hsug = h; ARIADNE_TESTING_NOTE(hsug);
            ValidatedVectorMultivariateFunctionPatch flow(integrator_ptr->flow_step(f,domx,suggest(hsug))); ARIADNE_TESTING_NOTE(flow);
            auto hused = flow.domain()[1].upper_bound(); ARIADNE_TESTING_NOTE(hused);
            ARIADNE_TESTING_EXPECT_LESS(flow.domain()[1].upper_bound(), h);
            if (hsug==hused) { ARIADNE_TESTING_NOTE("Integrator::flow_step(f,dom,suggest("<<h<<") used step-size "<<hused<<"; suggested step-size changed to actual step-size."); }
            else if (hsug==h) { ARIADNE_TESTING_NOTE("Integrator::flow_step(f,dom,suggest("<<h<<") used step-size "<<hused<<"; suggested step-size kept at requested step-size."); }
            else { { ARIADNE_TESTING_NOTE("Integrator::flow_step(f,d,suggest("<<h<<") used step-size "<<hused<<"; suggested step-size changed to "<<hsug<<", which is neither equal to requested step-size nor actual step-size."); ARIADNE_TESTING_PROVE(false); }; }
        } catch (IncompleteFlowException const& e) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,suggest("<<h<<")) caught IncompleteFlowException."; this->ctx().warning(message.str()); } while(false);
            ARIADNE_TESTING_EXPECT_LESS(e.computed_model().domain()[1].upper_bound(), h);
        } catch (FlowBoundsException const&) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,suggest("<<h<<")) threw FlowBoundsException, when a smaller step-size should be tried."; this->ctx().warning(message.str()); } while(false);
        } catch (FlowTimeStepException const&) {
            do { std::ostringstream message; message << "Integrator::flow_step(f,d,suggest("<<h<<")) threw FlowTimeStepException, when a smaller step-size should be tried."; this->ctx().warning(message.str()); } while(false);
        } catch (DivideByZeroException const&) {
            { ARIADNE_TESTING_NOTE("Integrator::flow_step(f,d,suggest("<<h<<")) threw DivideByZeroException, when a smaller step-size should be tried."); ARIADNE_TESTING_PROVE(false); };
        }

    }

};
class TestIntegratorExtras : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;


Void test_taylor_series_step_maximum_error() {
    auto x=ValidatedScalarMultivariateFunction::coordinates(2);
    ValidatedVectorMultivariateFunction f({
        x[1]+pow(x[0],3u)-pow(x[0],5u)+pow(x[0],3u)*pow(x[1],4u),
        -x[0]+pow(x[1],3u)+pow(x[0],4u)*pow(x[1],3u)-pow(x[1],7u)
    });
    ExactBoxType domain={
        {-0.25_x,0.25_x},
        {-0.25_x,0.25_x}
    };
    StepSizeType h=0.03125_x;

    TaylorSeriesIntegrator strict_integrator(step_maximum_error=1e-8,order=8);
    ARIADNE_TESTING_PROVE(strict_integrator.has_step_maximum_error());
    ARIADNE_TESTING_EXPECT_THROWS(FlowTimeStepException, strict_integrator.flow_step(f,domain,h));

    TaylorSeriesIntegrator loose_integrator(step_maximum_error=1e-2,order=8);
    FlowStepModelType flow(loose_integrator.flow_step(f,domain,h)); ARIADNE_TESTING_NOTE(flow);
    ARIADNE_TESTING_EXPECT_LESS(flow.error(), 1e-2_pr);
}

Void test_bounder_time_dependent_flow() {
    auto t=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({t});
    ExactBoxType domain={{0.0_x,1.0_x}};
    StepSizeType h(0.125_x);
    TaylorSeriesBounderIntegrator integrator(step_maximum_error=1e-6,order=4);
    auto flow=integrator.flow_step(field,domain,suggest(h));
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(flow.domain()[1u].upper_bound(), FloatDP(h,dp));
    auto time=FloatDPBounds(flow.domain()[1u].upper_bound());
    Vector<FloatDPBounds> point({FloatDPBounds(0.5_x,dp),time});
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.5_x,dp)+sqr(time)/2,flow(point)[0u]));
}


};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_integrator");
    ThresholdSweeper<FloatDP> sweeper(DoublePrecision(),1e-10);
    auto integrator_v0 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "TaylorPicardIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, TaylorPicardIntegrator(step_maximum_error=1e-8, sweeper, lipschitz_tolerance=0.5_x, minimum_temporal_order=0, maximum_temporal_order=16)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_spiral);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_logistic);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_time_variant);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_time_variant_with_parameters);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v0, test_step_size);
    auto integrator_v1 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "GradedTaylorPicardIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, GradedTaylorPicardIntegrator(step_maximum_error=1e-8,order=8)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_spiral);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_logistic);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_time_variant);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_time_variant_with_parameters);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v1, test_step_size);
    auto integrator_v2 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "TaylorSeriesIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, TaylorSeriesIntegrator(sweeper, lipschitz_tolerance=0.5_x, order=6)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_spiral);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_logistic);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_time_variant);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_time_variant_with_parameters);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v2, test_step_size);
    auto integrator_v3 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "TaylorSeriesBounderIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, TaylorSeriesBounderIntegrator(step_maximum_error=1e-6,sweeper,lipschitz_tolerance=0.5_x,order=6)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_spiral);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_logistic);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_time_variant);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_time_variant_with_parameters);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v3, test_step_size);
    auto integrator_v4 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "GradedTaylorSeriesIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, GradedTaylorSeriesIntegrator(step_maximum_error=1e-6, sweeper, lipschitz_tolerance=0.5_x, minimum_spacial_order=1, minimum_temporal_order=4, maximum_spacial_order=4, maximum_temporal_order=8)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_spiral);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_logistic);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_time_variant);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_time_variant_with_parameters);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v4, test_step_size);
    auto integrator_v5 = Ariadne::Testing::fixture_variant<TestIntegrator>(ariadne_test_suite, "AffineIntegrator", [&](Ariadne::Testing::Context& ctx) { return TestIntegrator(ctx, AffineIntegrator(1, 6)); });
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v5, test_constant_derivative);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v5, test_quadratic_flow);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v5, test_linear);
    ARIADNE_TESTING_ADD_VARIANT_TEST(integrator_v5, test_spiral);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorExtras, test_taylor_series_step_maximum_error);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorExtras, test_bounder_time_dependent_flow);
    // Preserve the legacy logging-verbosity CTest invocation separately from
    // Ariadne Testing's own --full/--concise reporting flags.
    if (argc == 3 && std::string_view(argv[1]) == "--verbosity" &&
        std::string_view(argv[2]) == "12") {
        Logger::instance().configuration().set_verbosity(12u);
        return ARIADNE_TESTING_RUN_TESTS();
    }
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
