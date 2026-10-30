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

#include <memory>
#include <sstream>
#include <limits>
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
Void check_constant_flow(IntegratorInterface const& integrator) {
    auto one=EffectiveScalarMultivariateFunction::constant(1u,1);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({one});
    ExactBoxType domain={{-0.125_x,0.125_x}};
    StepSizeType h(0.03125_x);
    auto flow=integrator.flow_step(field,domain,suggest(h));
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto t=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    EffectiveVectorMultivariateFunction expected={x+t};
    ARIADNE_TEST_COMPARE(norm(flow-expected),<,1e-6_pr)
    std::ostringstream text;
    text << integrator;
    ARIADNE_TEST_ASSERT(not text.str().empty())
}

Void test_picard_configuration() {
    TaylorPicardIntegrator integrator(StepMaximumError(1e-6));
    integrator.set_minimum_temporal_order(2u);
    integrator.set_maximum_temporal_order(10u);
    integrator.set_step_maximum_error(ApproximateDouble(1e-7));
    ARIADNE_TEST_EQUALS(integrator.minimum_temporal_order(),2u)
    ARIADNE_TEST_EQUALS(integrator.maximum_temporal_order(),10u)
    ARIADNE_TEST_EQUALS(integrator.step_maximum_error(),cast_exact(ApproximateDouble(1e-7)))
    Nat too_large=static_cast<Nat>(std::numeric_limits<DegreeType>::max())+1u;
    ARIADNE_TEST_THROWS(integrator.set_minimum_temporal_order(too_large),std::runtime_error)
    ARIADNE_TEST_THROWS(integrator.set_maximum_temporal_order(too_large),std::runtime_error)
    ARIADNE_TEST_EQUALS(integrator.minimum_temporal_order(),2u)
    ARIADNE_TEST_EQUALS(integrator.maximum_temporal_order(),10u)
    std::unique_ptr<IntegratorInterface> copy(integrator.clone());
    check_constant_flow(*copy);
    ARIADNE_TEST_THROWS(TaylorPicardIntegrator(step_maximum_error=1e-6,
        ThresholdSweeper<FloatDP>(dp,1e-10),lipschitz_tolerance=0.0_x,
        minimum_temporal_order=2,maximum_temporal_order=8),std::runtime_error)
}

Void test_graded_picard_configuration() {
    GradedTaylorPicardIntegrator integrator(step_maximum_error=1e-6,order=6);
    integrator.set_order(8u);
    integrator.set_step_maximum_error(ApproximateDouble(1e-7));
    integrator.set_error_refinement_minimum_improvement_percentage(ApproximateDouble(25.0));
    integrator.set_maximum_error_refinement_iterations(2u);
    integrator.set_diagnostics(true);
    ARIADNE_TEST_EQUALS(integrator.order(),8u)
    ARIADNE_TEST_EQUALS(integrator.step_maximum_error(),cast_exact(ApproximateDouble(1e-7)))
    ARIADNE_TEST_EQUALS(integrator.error_refinement_minimum_improvement_percentage(),25.0_x)
    ARIADNE_TEST_EQUALS(integrator.maximum_error_refinement_iterations(),2u)
    Nat too_large=static_cast<Nat>(std::numeric_limits<DegreeType>::max())+1u;
    ARIADNE_TEST_THROWS(integrator.set_order(too_large),std::runtime_error)
    ARIADNE_TEST_EQUALS(integrator.order(),8u)
    check_constant_flow(integrator);
    GradedTaylorPicardIntegrator swept(StepMaximumError(ApproximateDouble(1e-6)),
        Order(6u),StepSweepThreshold(ApproximateDouble(1e-12)));
    check_constant_flow(swept);
}

Void test_series_configuration() {
    ThresholdSweeper<FloatDP> sweeper(dp,1e-10);
    TaylorSeriesIntegrator integrator(sweeper,lipschitz_tolerance=0.5_x,order=4);
    ARIADNE_TEST_ASSERT(not integrator.has_step_maximum_error())
    ARIADNE_TEST_THROWS(integrator.step_maximum_error(),std::runtime_error)
    integrator.set_order(6u);
    integrator.set_sweeper(ThresholdSweeper<FloatDP>(dp,1e-12));
    check_constant_flow(integrator);
    integrator.set_step_maximum_error(ApproximateDouble(1e-7));
    ARIADNE_TEST_ASSERT(integrator.has_step_maximum_error())
    ARIADNE_TEST_EQUALS(integrator.order(),6u)
    ARIADNE_TEST_EQUALS(integrator.step_maximum_error(),cast_exact(ApproximateDouble(1e-7)))
    check_constant_flow(integrator);
    TaylorSeriesBounderIntegrator bounder_integrator(step_maximum_error=1e-6,order=6);
    check_constant_flow(bounder_integrator);
}

Void test_graded_series_constructors_and_configuration() {
    ThresholdSweeper<FloatDP> sweeper(dp,1e-10);
    GradedTaylorSeriesIntegrator defaults(step_maximum_error=1e-6);
    GradedTaylorSeriesIntegrator with_sweeper(step_maximum_error=1e-6,sweeper);
    GradedTaylorSeriesIntegrator with_lipschitz(step_maximum_error=1e-6,sweeper,lipschitz_tolerance=0.5_x);
    GradedTaylorSeriesIntegrator with_order(step_maximum_error=1e-6,sweeper,lipschitz_tolerance=0.5_x,maximum_temporal_order=8);
    with_order.set_minimum_spacial_order(1u);
    with_order.set_maximum_spacial_order(3u);
    with_order.set_minimum_temporal_order(2u);
    with_order.set_maximum_temporal_order(8u);
    with_order.set_sweeper(ThresholdSweeper<FloatDP>(dp,1e-12));
    with_order.set_step_maximum_error(ApproximateDouble(1e-7));
    ARIADNE_TEST_EQUALS(with_order.minimum_spacial_order(),1u)
    ARIADNE_TEST_EQUALS(with_order.maximum_spacial_order(),3u)
    ARIADNE_TEST_EQUALS(with_order.minimum_temporal_order(),2u)
    ARIADNE_TEST_EQUALS(with_order.maximum_temporal_order(),8u)
    ARIADNE_TEST_EQUALS(with_order.step_maximum_error(),cast_exact(ApproximateDouble(1e-7)))
    check_constant_flow(defaults);
    check_constant_flow(with_sweeper);
    check_constant_flow(with_lipschitz);
    check_constant_flow(with_order);
}

Void test_common_bounder_and_factory_interfaces() {
    TaylorSeriesIntegrator integrator(step_maximum_error=1e-6,order=6);
    EulerBounder bounder;
    integrator.set_bounder(bounder);
    std::ostringstream text;
    text << integrator.bounder();
    ARIADNE_TEST_EQUALS(text.str(),"EulerBounder")
    ValidatedFunctionPatchFactory factory(make_taylor_function_patch_factory(ThresholdSweeper<FloatDP>(dp,1e-12)));
    integrator.set_function_factory(factory);
    ExactBoxType domain={{-0.125_x,0.125_x}};
    ExactBoxType parameters={{-0.125_x,0.125_x}};
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto a=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction parameterised=EffectiveVectorMultivariateFunction({x+a});
    auto xp=EffectiveScalarMultivariateFunction::coordinate(3u,0u);
    auto t=EffectiveScalarMultivariateFunction::coordinate(3u,1u);
    auto p=EffectiveScalarMultivariateFunction::coordinate(3u,2u);
    ValidatedVectorMultivariateFunction timed=EffectiveVectorMultivariateFunction({xp+t+p});
    StepSizeType h(0.03125_x);
    auto first=integrator.flow_bounds(parameterised,domain,parameters,suggest(h));
    auto second=integrator.flow_bounds(timed,domain,StepSizeType(0),parameters,suggest(h));
    for(auto const& result_bounds:{first,second}) {
        ARIADNE_TEST_COMPARE(result_bounds.first,<=,h)
        ARIADNE_TEST_ASSERT(definitely(is_bounded(result_bounds.second)))
    }
    auto zero=EffectiveScalarMultivariateFunction::constant(1u,0);
    ValidatedVectorMultivariateFunction stationary=EffectiveVectorMultivariateFunction({zero});
    IntegratorInterface const& interface=integrator;
    auto automatic=interface.flow_step(stationary,domain);
    ARIADNE_TEST_EQUALS(automatic.argument_size(),2u)
    ARIADNE_TEST_EQUALS(automatic.domain()[1u].upper_bound(),1.0_x)
    check_constant_flow(integrator);
}

Void test_partial_flow_payload_and_flow_tube_printing() {
    TaylorPicardIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain={{-0.125_x,0.125_x}};
    auto zero=EffectiveScalarMultivariateFunction::constant(1u,0);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({zero});
    StepSizeType h(0.03125_x);
    auto flow=integrator.flow_step(field,domain,suggest(h));
    FlowBoundsException unavailable("bounds unavailable");
    ARIADNE_TEST_EQUALS(std::string(unavailable.what()),"bounds unavailable")
    IncompleteFlowException partial("partial flow",flow);
    IncompleteFlowException copied(partial);
    ARIADNE_TEST_EQUALS(std::string(copied.what()),"partial flow")
    ARIADNE_TEST_EQUALS(copied.computed_model().argument_size(),flow.argument_size())
    ARIADNE_TEST_EQUALS(copied.computed_model().result_size(),flow.result_size())
    ARIADNE_TEST_COMPARE(norm(copied.computed_model()-flow),<,1e-6_pr)
    FlowModelType tube;
    std::ostringstream empty;
    empty << tube;
    ARIADNE_TEST_EQUALS(empty.str(),"[  ]")
    tube.append(flow);
    tube.append(flow);
    std::ostringstream two_steps;
    two_steps << tube;
    ARIADNE_TEST_ASSERT(two_steps.str().find(",\n")!=std::string::npos)
}
}

int main() {
    ARIADNE_TEST_CALL(test_picard_configuration())
    ARIADNE_TEST_CALL(test_graded_picard_configuration())
    ARIADNE_TEST_CALL(test_series_configuration())
    ARIADNE_TEST_CALL(test_graded_series_constructors_and_configuration())
    ARIADNE_TEST_CALL(test_common_bounder_and_factory_interfaces())
    ARIADNE_TEST_CALL(test_partial_flow_payload_and_flow_tube_printing())
    return ARIADNE_TEST_FAILURES;
}
