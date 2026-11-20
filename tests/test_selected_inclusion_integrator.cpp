/***************************************************************************
 *            test_selected_inclusion_integrator.cpp
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

/** @file tests/test_selected_inclusion_integrator.cpp
 *  @brief Verify fixed and independent reach/evolve input approximations.
 */

#include <stdexcept>

#include "function/formula.tpl.hpp"
#include "function/taylor_function.hpp"
#include "integration/integrator.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestSelectedInclusionIntegrator : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    EffectiveVectorMultivariateFunction field() {
        auto x = EffectiveFormula::coordinate(0u);
        auto u = EffectiveFormula::coordinate(1u);
        return EffectiveVectorMultivariateFormulaFunction(2u, Vector<EffectiveFormula>{x + u});
    }

    BoxDomainType inputs() { return BoxDomainType({{-0.125_x, 0.125_x}}); }

    Void test_valid_and_invalid_strategy_pairs() {
        TaylorSeriesIntegrator base(step_maximum_error = 1e-6, order = 8);
        InclusionIntegratorFactory factory(base);
        auto f = field();
        auto U = inputs();
        for (const auto& reach: {InputApproximation(ZeroApproximation()),
                                 InputApproximation(ConstantApproximation())}) {
            for (const auto& evolve: {InputApproximation(ZeroApproximation()),
                                      InputApproximation(ConstantApproximation()),
                                      InputApproximation(AffineApproximation()),
                                      InputApproximation(SinusoidalApproximation()),
                                      InputApproximation(PiecewiseApproximation())}) {
                SelectedInclusionIntegrator selected(factory, f, U, reach, evolve);
                ARIADNE_TESTING_EXPECT_EQUAL(selected.reach_index(), reach.reference().index());
                ARIADNE_TESTING_EXPECT_EQUAL(selected.evolve_index(), evolve.reference().index());
            }
        }
        for (const auto& invalid_reach: {InputApproximation(AffineApproximation()),
                                         InputApproximation(SinusoidalApproximation()),
                                         InputApproximation(PiecewiseApproximation())}) {
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                SelectedInclusionIntegrator(factory, f, U, invalid_reach, InputApproximation(ZeroApproximation())));
        }
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            SelectedInclusionIntegrator(factory, f, BoxDomainType(0u),
                                        InputApproximation(ZeroApproximation()),
                                        InputApproximation(ZeroApproximation())));
    }

    Void test_distinct_approximations_and_endpoint() {
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        InclusionIntegratorFactory factory(base);
        SelectedInclusionIntegrator integrator(factory, field(), inputs(),
                                               InputApproximation(ZeroApproximation()),
                                               InputApproximation(AffineApproximation()));
        BoxDomainType domain({{0.0_x, 0.125_x}});
        UpperBoxType bound({{-1.0_x, 1.0_x}});
        auto initial = base.function_factory().create_identity(domain);
        StepSizeType h(0.03125_x);
        auto step = integrator.step(domain, initial, bound, TimeStepType(0u), h);

        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models[0u].result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models[0u].argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.evolve_model.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.evolve_model.argument_size(), 3u);
        auto expected = FloatDPBounds(0.0625_x, dp) * exp(FloatDPBounds(h, dp));
        UpperIntervalType expected_interval(expected.lower(), expected.upper());
        ARIADNE_TESTING_PROVE(refines(expected_interval, step.evolve_model.range()[0u]));
        const auto boxes = step.bounding_step(h);
        ARIADNE_TESTING_EXPECT_EQUAL(boxes.time_step, h);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(boxes.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(boxes.evolve)));
    }

    Void test_same_method_and_nonrigorous_endpoint_strategy() {
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        InclusionIntegratorFactory factory(base);
        BoxDomainType domain({{0.0_x, 0.125_x}});
        UpperBoxType bound({{-1.0_x, 1.0_x}});
        auto initial = base.function_factory().create_identity(domain);
        StepSizeType h(0.03125_x);
        for (const auto& approximation: {InputApproximation(ZeroApproximation()),
                                         InputApproximation(ConstantApproximation())}) {
            SelectedInclusionIntegrator selected(factory, field(), inputs(), approximation, approximation);
            auto step = selected.step(domain, initial, bound, TimeStepType(0u), h);
            ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models.size(), 1u);
            ARIADNE_TESTING_PROVE(definitely(is_bounded(step.bounding_step(h).reach)));
            ARIADNE_TESTING_PROVE(definitely(is_bounded(step.bounding_step(h).evolve)));
        }
        for (const auto& approximation: {InputApproximation(SinusoidalApproximation()),
                                         InputApproximation(PiecewiseApproximation())}) {
            SelectedInclusionIntegrator selected(factory, field(), inputs(),
                                                 InputApproximation(ZeroApproximation()), approximation);
            auto step = selected.step(domain, initial, bound, TimeStepType(0u), h);
            ARIADNE_TESTING_PROVE(not approximation.is_reach_rigorous());
            ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models.size(), 1u);
            ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolve_model.range())));
        }
    }

    Void test_invalid_time_step() {
        TaylorSeriesIntegrator base(step_maximum_error = 1e-6, order = 8);
        InclusionIntegratorFactory factory(base);
        SelectedInclusionIntegrator selected(factory, field(), inputs(),
                                             InputApproximation(ZeroApproximation()),
                                             InputApproximation(ConstantApproximation()));
        BoxDomainType domain({{0.0_x, 0.125_x}});
        auto initial = base.function_factory().create_identity(domain);
        UpperBoxType bound({{-1.0_x, 1.0_x}});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            selected.step(domain, initial, bound, TimeStepType(0u), StepSizeType(0u)));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_selected_inclusion_integrator");
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_valid_and_invalid_strategy_pairs);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_distinct_approximations_and_endpoint);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_same_method_and_nonrigorous_endpoint_strategy);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_invalid_time_step);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
