/***************************************************************************
 *            test_preconditioned_integrator.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file tests/test_preconditioned_integrator.cpp
 *  @brief Test preconditioned integrator behavior.
 */

#include <cmath>
#include <memory>
#include <sstream>
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "logging/logging.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {
class TestPreconditionedIntegrator : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

Void check_factorisation(PreconditionedTaylorSeriesState const& state,
                         ValidatedVectorMultivariateFunctionPatch const& original,
                         Vector<FloatDPBounds> const& point) {
    auto expected=original(point);
    auto local=state.normalised_mapping()(point);
    ARIADNE_TESTING_EXPECT_EQUAL(state.parameter_domain().dimension(),original.argument_size());
    ARIADNE_TESTING_EXPECT_EQUAL(state.local_domain().dimension(),original.result_size());
    for(SizeType i=0u;i<expected.size();++i) {
        FloatDPBounds reconstructed(state.centre()[i]);
        for(SizeType j=0u;j<local.size();++j) {
            reconstructed+=FloatDPBounds(state.linear_map()[i][j])*local[j];
        }
        ARIADNE_TESTING_PROVE(refines(expected[i],reconstructed));
    }
}

Void test_preconditioning_geometry() {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain={{-1.0_x,1.0_x},{-1.0_x,1.0_x}};
    auto const& factory=integrator.function_factory();
    auto id=factory.create_identity(domain);
    auto mapping=factory.create_zeros(2u,domain);
    mapping[0u]=FloatDPBounds(1,dp)+id[0u]+FloatDPBounds(3,dp)*id[1u];
    mapping[1u]=FloatDPBounds(2,dp)+FloatDPBounds(2,dp)*id[0u]+hlf(id[1u]);
    for(auto mode:{TaylorSeriesPreconditioning::IDENTITY,TaylorSeriesPreconditioning::QR}) {
        integrator.set_preconditioning(mode);
        auto state=integrator.precondition(mapping);
        check_factorisation(state,mapping,Vector<FloatDPBounds>({1,-1},dp));
        check_factorisation(state,mapping,Vector<FloatDPBounds>({0,0},dp));
        for(SizeType i=0u;i<2u;++i) {
            for(SizeType j=0u;j<2u;++j) {
                double dot=0.0;
                for(SizeType k=0u;k<2u;++k) {
                    dot+=state.linear_map()[k][i].get_d()*state.linear_map()[k][j].get_d();
                }
                ARIADNE_TESTING_EXPECT_WITHIN(dot,(i==j ? 1.0 : 0.0),1e-12);
            }
        }
        std::ostringstream text;
        text << static_cast<IntegratorInterface const&>(integrator);
        ARIADNE_TESTING_PROVE(text.str().find(mode==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY")!=std::string::npos);
        std::unique_ptr<PreconditionedGradedTaylorSeriesIntegrator> copy(integrator.clone());
        ARIADNE_TESTING_PROVE(copy->preconditioning()==mode);
    }
    // Extra enclosure parameters intentionally retain identity orientation.
    ExactBoxType parameter_domain(3u,ExactIntervalType(-1.0_x,1.0_x));
    auto parameters=factory.create_identity(parameter_domain);
    auto rectangular=factory.create_zeros(2u,parameter_domain);
    rectangular[0u]=parameters[0u]+parameters[2u];
    rectangular[1u]=parameters[1u];
    check_factorisation(integrator.precondition(rectangular),rectangular,
                        Vector<FloatDPBounds>({1,-1,1},dp));
}

Void test_rank_deficient_geometry() {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain(2u,ExactIntervalType(-1.0_x,1.0_x));
    auto const& factory=integrator.function_factory();
    auto constant=factory.create_zeros(2u,domain);
    constant[0u]=factory.create_constant(domain,1);
    constant[1u]=factory.create_constant(domain,2);
    auto state=integrator.precondition(constant);
    check_factorisation(state,constant,Vector<FloatDPBounds>({0,0},dp));
    ARIADNE_TESTING_EXPECT_EQUAL(state.linear_map()[0u][0u],FloatDP(1,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(state.linear_map()[1u][1u],FloatDP(1,dp));
    // The coefficients are finite, but forming a squared norm can overflow.
    // The preconditioner must retain a usable identity orientation.
    ExactBoxType large_domain={{-1.0_x,1.0_x}};
    auto large=factory.create_identity(large_domain);
    large[0u]=large[0u]*FloatDPBounds(0x1p1020_x,dp);
    auto large_state=integrator.precondition(large);
    ARIADNE_TESTING_EXPECT_EQUAL(large_state.linear_map()[0u][0u],FloatDP(1,dp));
    check_factorisation(large_state,large,Vector<FloatDPBounds>({FloatDPBounds(1,dp)}));
}

Void check_rotation_endpoint(PreconditionedTaylorSeriesStep const& step,StepSizeType const& elapsed) {
    Vector<FloatDPBounds> point({1,1},dp);
    auto value=step.evolved_mapping()(point);
    FloatDPBounds t(elapsed,dp);
    Vector<FloatDPBounds> expected({cos(t)-sin(t),sin(t)+cos(t)});
    for(SizeType i=0u;i<2u;++i) {
        ARIADNE_TESTING_PROVE(refines(expected[i],value[i]));
    }
    ARIADNE_TESTING_EXPECT_EQUAL(step.final_state().centre().size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(step.final_state().parameter_domain().dimension(),2u);
    ARIADNE_TESTING_PROVE(definitely(is_bounded(step.flowpipe_mapping().range())));
}

Void test_rotating_flow_and_carried_state() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({-y,x});
    ExactBoxType domain(2u,ExactIntervalType(0.875_x,1.125_x));
    StepSizeType h(0.03125_x);
    for(auto mode:{TaylorSeriesPreconditioning::IDENTITY,TaylorSeriesPreconditioning::QR}) {
      for(Bool diagnostics:{false,true}) {
        ARIADNE_TESTING_NOTE("Two rotation steps, mode " << (mode==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY") << ", diagnostics " << diagnostics);
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
        integrator.set_preconditioning(mode);
        integrator.set_diagnostics(diagnostics);
        integrator.set_carried_expansion_diagnostics(diagnostics);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.diagnostics(),diagnostics);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.carried_expansion_diagnostics(),diagnostics);
        auto initial=integrator.function_factory().create_identity(domain);
        auto first=integrator.step(field,integrator.precondition(initial),suggest(h));
        ARIADNE_TESTING_EXPECT_GREATER(first.time_step(), StepSizeType(0));
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(first.time_step(), h);
        check_rotation_endpoint(first,first.time_step());
        auto second=integrator.step(field,first.final_state(),suggest(h));
        check_rotation_endpoint(second,first.time_step()+second.time_step());
      }
    }
}

Void test_physical_flow_step() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({x});
    UpperBoxType bound={{0.5_x,2.0_x}};
    StepSizeType h(0.03125_x);
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    // The point domain exercises the zero-radius fallback; the interval
    // exercises normalisation and conversion back to physical coordinates.
    for(auto domain:{ExactBoxType{{0.875_x,1.125_x}},ExactBoxType{{1.0_x,1.0_x}}}) {
        auto flow=integrator.flow_step(field,domain,h,bound);
        ARIADNE_TESTING_EXPECT_EQUAL(flow.domain(),product(domain,ExactIntervalType(0,h)));
        for(auto time:{StepSizeType(0),h}) {
            Vector<FloatDPBounds> point({FloatDPBounds(1,dp),FloatDPBounds(time,dp)});
            ARIADNE_TESTING_PROVE(refines(exp(FloatDPBounds(time,dp)),flow(point)[0u]));
        }
    }
    ExactBoxType domain={{0.875_x,1.125_x}};
    integrator.set_minimum_spacial_order(1u);
    integrator.set_maximum_spacial_order(1u);
    integrator.set_minimum_temporal_order(1u);
    integrator.set_maximum_temporal_order(1u);
    integrator.set_step_maximum_error(1e-16);
    ARIADNE_TESTING_EXPECT_THROWS(FlowTimeStepException, integrator.flow_step(field,domain,h,bound));
    ValidatedVectorMultivariateFunction wrong_result=EffectiveVectorMultivariateFunction({x,x});
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.flow_step(wrong_result,domain,h,bound));
    auto extra=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction wrong_arguments=EffectiveVectorMultivariateFunction({extra});
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.flow_step(wrong_arguments,domain,h,bound));
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.flow_step(field,domain,h,UpperBoxType(0u)));
}

Void test_step_dimension_contracts() {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain={{-1.0_x,1.0_x}};
    auto mapping=integrator.function_factory().create_identity(domain);
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({x});
    ValidatedVectorMultivariateFunction wrong_result=EffectiveVectorMultivariateFunction({x,x});
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction wrong_arguments=EffectiveVectorMultivariateFunction({y});
    auto state=integrator.precondition(mapping);
    auto h=suggest(StepSizeType(0.03125_x));
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(wrong_result,state,h));
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(wrong_arguments,state,h));
    Vector<FloatDP> centre({FloatDP(0,dp)});
    PreconditionedTaylorSeriesState wrong_rows(centre,Matrix<FloatDP>(0u,1u,dp),domain,mapping);
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(field,wrong_rows,h));
    PreconditionedTaylorSeriesState wrong_columns(centre,Matrix<FloatDP>(1u,0u,dp),domain,mapping);
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(field,wrong_columns,h));
    PreconditionedTaylorSeriesState wrong_domain(centre,Matrix<FloatDP>::identity(1u,dp),ExactBoxType(),mapping);
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(field,wrong_domain,h));
    auto wrong_mapping=integrator.function_factory().create_zeros(2u,domain);
    PreconditionedTaylorSeriesState wrong_image(centre,Matrix<FloatDP>::identity(1u,dp),domain,wrong_mapping);
    ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step(field,wrong_image,h));
}

Void test_stationary_zero_order_step() {
    auto zero=EffectiveScalarMultivariateFunction::constant(1u,0);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({zero});
    for(Bool diagnostics:{false,true}) {
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
        integrator.set_minimum_spacial_order(0u);
        integrator.set_maximum_spacial_order(0u);
        integrator.set_minimum_temporal_order(0u);
        integrator.set_maximum_temporal_order(0u);
        integrator.set_diagnostics(diagnostics);
        ExactBoxType domain={{-1.0_x,1.0_x}};
        auto initial=integrator.function_factory().create_constants(domain,Vector<FloatDPBounds>({FloatDPBounds(1,dp)}));
        auto h=StepSizeType(0.03125_x);
        auto step=integrator.step(field,integrator.precondition(initial),suggest(h));
        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step(),h);
        auto value=step.evolved_mapping()(Vector<FloatDPBounds>({FloatDPBounds(0,dp)}));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(1,dp),value[0u]));
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(step.evolved_mapping().error(), 1e-6_pr);
    }
}

Void test_adaptive_orders_and_step_rejection() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({x});
    ExactBoxType domain={{-1.0_x,1.0_x}};
    StepSizeType h(0.125_x);
    for(unsigned verbosity:{0u,12u}) {
      Logger::instance().configuration().set_verbosity(verbosity);
      for(unsigned configuration:{0u,1u,2u}) {
       for(Bool diagnostics:{false,true}) {
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
        integrator.set_diagnostics(diagnostics);
        integrator.set_minimum_spacial_order(1u);
        integrator.set_maximum_spacial_order(configuration==1u ? 2u : 1u);
        integrator.set_minimum_temporal_order(1u);
        integrator.set_maximum_temporal_order(configuration==2u ? 2u : 1u);
        auto initial=integrator.function_factory().create_identity(domain);
        auto step=integrator.step(field,integrator.precondition(initial),suggest(h));
        ARIADNE_TESTING_EXPECT_GREATER(step.time_step(), StepSizeType(0));
        ARIADNE_TESTING_EXPECT_LESS(step.time_step(), h);
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(step.evolved_mapping().error(), 1e-6_pr);
        auto value=step.evolved_mapping()(Vector<FloatDPBounds>({FloatDPBounds(0.5_x,dp)}));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.5_x,dp)*exp(FloatDPBounds(step.time_step(),dp)),value[0u]));
       }
      }
    }
}
};

}


int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_preconditioned_integrator");
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_preconditioning_geometry);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_rank_deficient_geometry);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_rotating_flow_and_carried_state);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_physical_flow_step);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_step_dimension_contracts);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_stationary_zero_order_step);
    ARIADNE_TESTING_ADD_TEST(TestPreconditionedIntegrator, test_adaptive_orders_and_step_rejection);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
