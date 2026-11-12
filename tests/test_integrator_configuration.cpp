/***************************************************************************
 *            test_integrator_configuration.cpp
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

/** @file tests/test_integrator_configuration.cpp
 *  @brief Test integrator configuration behavior.
 */

#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include "function/taylor_function.hpp"
#include "integration/integrator.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {
class TestIntegratorConfiguration : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void check_constant_flow(const IntegratorInterface& integrator) {
        auto one = EffectiveScalarMultivariateFunction::constant(1u, 1);
        ValidatedVectorMultivariateFunction field = EffectiveVectorMultivariateFunction({one});
        ExactBoxType domain = {
            {-0.125_x, 0.125_x}
        };
        StepSizeType h(0.03125_x);
        auto flow = integrator.flow_step(field, domain, suggest(h));
        auto x = EffectiveScalarMultivariateFunction::coordinate(2u, 0u);
        auto t = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        EffectiveVectorMultivariateFunction expected = {x + t};
        ARIADNE_TESTING_EXPECT_LESS(norm(flow - expected), 1e-6_pr);
        std::ostringstream text;
        text << integrator;
        ARIADNE_TESTING_PROVE(not text.str().empty());
    }

    Void test_picard_configuration() {
        TaylorPicardIntegrator integrator(StepMaximumError(1e-6));
        integrator.set_minimum_temporal_order(2u);
        integrator.set_maximum_temporal_order(10u);
        integrator.set_step_maximum_error(ApproximateDouble(1e-7));
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.minimum_temporal_order(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.maximum_temporal_order(), 10u);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.step_maximum_error(), cast_exact(ApproximateDouble(1e-7)));
        Nat too_large = static_cast<Nat>(std::numeric_limits<DegreeType>::max()) + 1u;
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.set_minimum_temporal_order(too_large));
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.set_maximum_temporal_order(too_large));
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.minimum_temporal_order(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.maximum_temporal_order(), 10u);
        std::unique_ptr<IntegratorInterface> copy(integrator.clone());
        check_constant_flow(*copy);
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, TaylorPicardIntegrator(step_maximum_error = 1e-6,
                                                                                 ThresholdSweeper<FloatDP>(dp, 1e-10), lipschitz_tolerance = 0.0_x,
                                                                                 minimum_temporal_order = 2, maximum_temporal_order = 8));
    }

    Void test_graded_picard_configuration() {
        GradedTaylorPicardIntegrator integrator(step_maximum_error = 1e-6, order = 6);
        integrator.set_order(8u);
        integrator.set_step_maximum_error(ApproximateDouble(1e-7));
        integrator.set_error_refinement_minimum_improvement_percentage(ApproximateDouble(25.0));
        integrator.set_maximum_error_refinement_iterations(2u);
        integrator.set_diagnostics(true);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.order(), 8u);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.step_maximum_error(), cast_exact(ApproximateDouble(1e-7)));
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.error_refinement_minimum_improvement_percentage(), 25.0_x);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.maximum_error_refinement_iterations(), 2u);
        Nat too_large = static_cast<Nat>(std::numeric_limits<DegreeType>::max()) + 1u;
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.set_order(too_large));
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.order(), 8u);
        check_constant_flow(integrator);
        GradedTaylorPicardIntegrator swept(StepMaximumError(ApproximateDouble(1e-6)),
                                           Order(6u), StepSweepThreshold(ApproximateDouble(1e-12)));
        check_constant_flow(swept);
    }

    Void test_series_configuration() {
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        TaylorSeriesIntegrator integrator(sweeper, lipschitz_tolerance = 0.5_x, order = 4);
        ARIADNE_TESTING_PROVE(not integrator.has_step_maximum_error());
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator.step_maximum_error());
        integrator.set_order(6u);
        integrator.set_sweeper(ThresholdSweeper<FloatDP>(dp, 1e-12));
        check_constant_flow(integrator);
        integrator.set_step_maximum_error(ApproximateDouble(1e-7));
        ARIADNE_TESTING_PROVE(integrator.has_step_maximum_error());
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.order(), 6u);
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.step_maximum_error(), cast_exact(ApproximateDouble(1e-7)));
        check_constant_flow(integrator);
        TaylorSeriesBounderIntegrator bounder_integrator(step_maximum_error = 1e-6, order = 6);
        check_constant_flow(bounder_integrator);
    }

    Void test_graded_series_constructors_and_configuration() {
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        GradedTaylorSeriesIntegrator defaults(step_maximum_error = 1e-6);
        GradedTaylorSeriesIntegrator with_sweeper(step_maximum_error = 1e-6, sweeper);
        GradedTaylorSeriesIntegrator with_lipschitz(step_maximum_error = 1e-6, sweeper, lipschitz_tolerance = 0.5_x);
        GradedTaylorSeriesIntegrator with_order(step_maximum_error = 1e-6, sweeper, lipschitz_tolerance = 0.5_x, maximum_temporal_order = 8);
        with_order.set_minimum_spacial_order(1u);
        with_order.set_maximum_spacial_order(3u);
        with_order.set_minimum_temporal_order(2u);
        with_order.set_maximum_temporal_order(8u);
        with_order.set_sweeper(ThresholdSweeper<FloatDP>(dp, 1e-12));
        with_order.set_step_maximum_error(ApproximateDouble(1e-7));
        ARIADNE_TESTING_EXPECT_EQUAL(with_order.minimum_spacial_order(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(with_order.maximum_spacial_order(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(with_order.minimum_temporal_order(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(with_order.maximum_temporal_order(), 8u);
        ARIADNE_TESTING_EXPECT_EQUAL(with_order.step_maximum_error(), cast_exact(ApproximateDouble(1e-7)));
        check_constant_flow(defaults);
        check_constant_flow(with_sweeper);
        check_constant_flow(with_lipschitz);
        check_constant_flow(with_order);
    }

    Void test_common_bounder_and_factory_interfaces() {
        TaylorSeriesIntegrator integrator(step_maximum_error = 1e-6, order = 6);
        EulerBounder bounder;
        integrator.set_bounder(bounder);
        std::ostringstream text;
        text << integrator.bounder();
        ARIADNE_TESTING_EXPECT_EQUAL(text.str(), "EulerBounder");
        ValidatedFunctionPatchFactory factory(make_taylor_function_patch_factory(ThresholdSweeper<FloatDP>(dp, 1e-12)));
        integrator.set_function_factory(factory);
        ExactBoxType domain = {
            {-0.125_x, 0.125_x}
        };
        ExactBoxType parameters = {
            {-0.125_x, 0.125_x}
        };
        auto x = EffectiveScalarMultivariateFunction::coordinate(2u, 0u);
        auto a = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        ValidatedVectorMultivariateFunction parameterised = EffectiveVectorMultivariateFunction({x + a});
        auto xp = EffectiveScalarMultivariateFunction::coordinate(3u, 0u);
        auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        auto p = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        ValidatedVectorMultivariateFunction timed = EffectiveVectorMultivariateFunction({xp + t + p});
        StepSizeType h(0.03125_x);
        auto first = integrator.flow_bounds(parameterised, domain, parameters, suggest(h));
        auto second = integrator.flow_bounds(timed, domain, StepSizeType(0), parameters, suggest(h));
        for (const auto& result_bounds: {first, second}) {
            ARIADNE_TESTING_EXPECT_LESS_EQUAL(result_bounds.first, h);
            ARIADNE_TESTING_PROVE(definitely(is_bounded(result_bounds.second)));
        }
        auto zero = EffectiveScalarMultivariateFunction::constant(1u, 0);
        ValidatedVectorMultivariateFunction stationary = EffectiveVectorMultivariateFunction({zero});
        const IntegratorInterface& interface = integrator;
        auto automatic = interface.flow_step(stationary, domain);
        ARIADNE_TESTING_EXPECT_EQUAL(automatic.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(automatic.domain()[1u].upper_bound(), 1.0_x);
        check_constant_flow(integrator);
    }

    Void test_partial_flow_payload_and_flow_tube_printing() {
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ExactBoxType domain = {
            {-0.125_x, 0.125_x}
        };
        auto zero = EffectiveScalarMultivariateFunction::constant(1u, 0);
        ValidatedVectorMultivariateFunction field = EffectiveVectorMultivariateFunction({zero});
        StepSizeType h(0.03125_x);
        auto flow = integrator.flow_step(field, domain, suggest(h));
        FlowBoundsException unavailable("bounds unavailable");
        ARIADNE_TESTING_EXPECT_EQUAL(std::string(unavailable.what()), "bounds unavailable");
        IncompleteFlowException partial("partial flow", flow);
        IncompleteFlowException copied(partial);
        ARIADNE_TESTING_EXPECT_EQUAL(std::string(copied.what()), "partial flow");
        ARIADNE_TESTING_EXPECT_EQUAL(copied.computed_model().argument_size(), flow.argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(copied.computed_model().result_size(), flow.result_size());
        ARIADNE_TESTING_EXPECT_LESS(norm(copied.computed_model() - flow), 1e-6_pr);
        FlowModelType tube;
        std::ostringstream empty;
        empty << tube;
        ARIADNE_TESTING_EXPECT_EQUAL(empty.str(), "[  ]");
        tube.append(flow);
        tube.append(flow);
        std::ostringstream two_steps;
        two_steps << tube;
        ARIADNE_TESTING_PROVE(two_steps.str().find(",\n") != std::string::npos);
    }

    Void test_flow_dimension_contracts() {
        auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        auto extra = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        ValidatedVectorMultivariateFunction field = EffectiveVectorMultivariateFunction({x});
        ValidatedVectorMultivariateFunction wrong_result = EffectiveVectorMultivariateFunction({x, x});
        ValidatedVectorMultivariateFunction wrong_arguments = EffectiveVectorMultivariateFunction({extra});
        ExactBoxType domain = {
            {-0.125_x, 0.125_x}
        };
        UpperBoxType bound = {
            {-1.0_x, 1.0_x}
        };
        StepSizeType h(0.03125_x);
        Interval<StepSizeType> time(0, h);
        std::vector<std::unique_ptr<IntegratorInterface>> integrators;
        integrators.push_back(std::make_unique<TaylorPicardIntegrator>(step_maximum_error = 1e-6));
        integrators.push_back(std::make_unique<GradedTaylorPicardIntegrator>(step_maximum_error = 1e-6, order = 4));
        integrators.push_back(std::make_unique<TaylorSeriesIntegrator>(step_maximum_error = 1e-6, order = 4));
        integrators.push_back(std::make_unique<GradedTaylorSeriesIntegrator>(step_maximum_error = 1e-6));
        integrators.push_back(std::make_unique<AffineIntegrator>(SpacialOrder(1), TemporalOrder(4)));
        for (const auto& integrator: integrators) {
            ARIADNE_TESTING_NOTE(*integrator << " dimension contracts");
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator->flow_step(wrong_result, domain, h, bound));
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator->flow_step(wrong_arguments, domain, h, bound));
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator->flow_step(wrong_result, domain, time, ExactBoxType(), bound));
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrator->flow_step(wrong_arguments, domain, time, ExactBoxType(), bound));
        }
        auto t = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        ValidatedVectorMultivariateFunction timed = EffectiveVectorMultivariateFunction({t});
        UpperBoxType wrong_bound(0u);
        for (SizeType i = 0u; i < integrators.size(); ++i) {
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrators[i]->flow_step(field, domain, h, wrong_bound));
            ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, integrators[i]->flow_step(timed, domain, time, ExactBoxType(), wrong_bound));
        }
        TaylorSeriesBounderIntegrator bounder(step_maximum_error = 1e-6, order = 4);
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, bounder.flow_step(wrong_result, domain, suggest(h)));
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, bounder.flow_step(wrong_arguments, domain, suggest(h)));
    }

    Void test_stationary_point_flow() {
        auto zero = EffectiveScalarMultivariateFunction::constant(1u, 0);
        ValidatedVectorMultivariateFunction field = EffectiveVectorMultivariateFunction({zero});
        TaylorSeriesIntegrator integrator(step_maximum_error = 1e-6, order = 2);
        ExactBoxType domain = {
            {0.0_x, 0.0_x}
        };
        UpperBoxType bound = {
            {0.0_x, 0.0_x}
        };
        StepSizeType h(0.125_x);
        auto flow = integrator.flow_step(field, domain, h, bound);
        auto value = flow(Vector<FloatDPBounds>({FloatDPBounds(0, dp), FloatDPBounds(h, dp)}));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0, dp), value[0u]));
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(flow.error(), 1e-6_pr);
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_integrator_configuration");
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_picard_configuration);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_graded_picard_configuration);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_series_configuration);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_graded_series_constructors_and_configuration);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_common_bounder_and_factory_interfaces);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_partial_flow_payload_and_flow_tube_printing);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_flow_dimension_contracts);
    ARIADNE_TESTING_ADD_TEST(TestIntegratorConfiguration, test_stationary_point_flow);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
