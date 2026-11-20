/***************************************************************************
 *            test_validated_integration_step.cpp
 *
 *  Copyright  2026  Luca Geretti
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Integration. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file test_validated_integration_step.cpp
 *  @brief Verify one numerical result contract for ordinary and set-valued flows.
 */

#include <stdexcept>

#include "function/symbolic_function.hpp"
#include "integration/integration_step_adapters.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestValidatedIntegrationStep : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_single_valued_flow_preserves_patch() {
        // x'=1, X0=[0,1/8], h=1/32. The entire tube and the
        // endpoint must be certified; the original time-dependent
        // flow patch must still be available after box projection.
        EffectiveVectorMultivariateFunction f = {
            EffectiveScalarMultivariateFunction::constant(1u, 1)};
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        TaylorSeriesIntegrator solver(step_maximum_error = 1e-6, order = 8);
        const auto flow = solver.flow_step(f, initial, h, bound);
        const auto step = validated_integration_step(
            flow, Interval<StepSizeType>(StepSizeType(0u), h));

        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0_x, 0.15625_x}}), cast_exact_box(step.reach()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.03125_x, 0.15625_x}}), cast_exact_box(step.evolve()))));
        ARIADNE_TESTING_PROVE(step.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models().size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models()[0u].domain(), flow.domain());
        ARIADNE_TESTING_PROVE(step.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(step.evolve_model().argument_size(), initial.dimension());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.03125_x, 0.15625_x}}),
            cast_exact_box(UpperBoxType(step.evolve_model().range())))));
        ARIADNE_TESTING_PROVE(step.has_representation());
        const auto* restored = step.representation_as<FlowStepModelType>();
        ARIADNE_TESTING_PROVE(restored != nullptr);
        ARIADNE_TESTING_EXPECT_EQUAL(restored->domain(), flow.domain());
        ARIADNE_TESTING_EXPECT_EQUAL(restored->argument_size(), 2u);
        ARIADNE_TESTING_PROVE(step.representation_as<SelectedInclusionStep>() == nullptr);

        // The box enclosures and the structured certificates must agree.
        // An undersized reach bound cannot be advertised as the hull of a
        // certified full-time flow; likewise for the endpoint map.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(
                h, UpperBoxType({{0.0625_x, 0.15625_x}}), step.evolve(),
                {}, step.reach_models(), step.evolve_model(), StepSizeType(0u), 1u));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(
                h, step.reach(), UpperBoxType({{0.0625_x, 0.15625_x}}),
                {}, step.reach_models(), step.evolve_model(), StepSizeType(0u), 1u));

        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            validated_integration_step(
                flow, Interval<StepSizeType>(StepSizeType(0u), StepSizeType(0.0625_x))));
        // The autonomous model's [0,h] time is local. Nonautonomous models
        // must not falsely reinterpret that time domain as absolute.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            validated_integration_step(
                flow, Interval<StepSizeType>(StepSizeType(0u), h),
                StepSizeType(0.5_x), std::nullopt, true));
    }

    Void test_piecewise_reach_sections_cover_every_time() {
        // Two certified restrictions of the SAME numerical flow provide a
        // multi-patch reach. They deliberately arrive in reverse temporal
        // order; the common contract must not depend on patch ordering.
        const EffectiveVectorMultivariateFunction field = {
            EffectiveScalarMultivariateFunction::constant(1u, 1)};
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        const StepSizeType half(0.015625_x);
        TaylorSeriesIntegrator solver(step_maximum_error = 1e-6, order = 8);
        const auto flow = solver.flow_step(field, initial, h, bound);
        const auto native = validated_integration_step(
            flow, Interval<StepSizeType>(StepSizeType(0u), h));

        const ValidatedVectorMultivariateFunctionPatch full = flow;
        const auto first_half = restriction(full, product(
            BoxDomainType(initial), IntervalDomainType(0.0_x, 0.015625_x)));
        const auto last_half = restriction(full, product(
            BoxDomainType(initial), IntervalDomainType(0.015625_x, 0.03125_x)));

        List<ValidatedVectorMultivariateFunctionPatch> pieces;
        pieces.append(last_half);
        pieces.append(first_half);
        const ValidatedIntegrationStep step(
            h, native.reach(), native.evolve(), {}, pieces,
            native.evolve_model(), StepSizeType(0u), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_sections(StepSizeType(0u)).size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_sections(half).size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_sections(h).size(), 1u);

        const auto left = step.reach_sections(StepSizeType(0.0078125_x));
        const auto right = step.reach_sections(StepSizeType(0.0234375_x));
        ARIADNE_TESTING_EXPECT_EQUAL(left.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(right.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(left[0u].domain(), BoxDomainType(initial));
        ARIADNE_TESTING_EXPECT_EQUAL(right[0u].domain(), BoxDomainType(initial));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0078125_x, 0.1328125_x}}),
            cast_exact_box(UpperBoxType(left[0u].range())))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0234375_x, 0.1484375_x}}),
            cast_exact_box(UpperBoxType(right[0u].range())))));

        // A partial tube must not pretend to provide certified time
        // sections for an entire step, even if its overall box is sound.
        List<ValidatedVectorMultivariateFunctionPatch> incomplete;
        incomplete.append(first_half);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(h, native.reach(), native.evolve(),
                                     {}, incomplete, native.evolve_model(),
                                     StepSizeType(0u), 1u));

        // The same local-time patch cannot be passed off as a reach
        // whose time coordinates are absolute around physical t=1/2.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(h, native.reach(), native.evolve(),
                                     {}, pieces, native.evolve_model(),
                                     StepSizeType(0.5_x), 1u));
    }

    Void test_set_valued_box_step_without_forced_structure() {
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        const BoxDomainType states({{-1.0_x, 1.0_x}});
        const BoxDomainType inputs({{-1.0_x, 1.0_x}});
        ValidatedInputFieldBoxAdapterDP field(f, inputs, states);
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        FirstOrderMultifunctionIntegrator solver;
        const auto native = solver.flow_step(field, initial, h, bound);
        const auto step = validated_integration_step(native);

        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}), cast_exact_box(step.reach()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}), cast_exact_box(step.evolve()))));
        ARIADNE_TESTING_PROVE(not step.has_reach_models());
        ARIADNE_TESTING_PROVE(not step.has_evolve_model());
        ARIADNE_TESTING_EXPECT_THROWS(std::logic_error, step.evolve_model());
        ARIADNE_TESTING_PROVE(not step.has_representation());
    }

    Void test_structured_inclusion_retains_parameter_maps() {
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        const BoxDomainType states({{-1.0_x, 1.0_x}});
        const BoxDomainType inputs({{-1.0_x, 1.0_x}});
        ValidatedInputFieldBoxAdapterDP field(f, inputs, states);
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator solver(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const auto native = solver.structured_step(field, initial, h, bound);
        const auto step = validated_integration_step(native);

        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}), cast_exact_box(step.reach()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}), cast_exact_box(step.evolve()))));
        ARIADNE_TESTING_PROVE(step.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(step.reach_models().size(), native.reach_models.size());
        ARIADNE_TESTING_PROVE(step.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(step.evolve_model().domain(), native.evolve_model.domain());
        ARIADNE_TESTING_EXPECT_EQUAL(
            step.evolve_model().argument_size(), native.evolve_model.argument_size());
        const auto* restored = step.representation_as<SelectedInclusionStep>();
        ARIADNE_TESTING_PROVE(restored != nullptr);
        ARIADNE_TESTING_EXPECT_EQUAL(restored->time_step, native.time_step);
        ARIADNE_TESTING_EXPECT_EQUAL(
            restored->evolve_model.domain(), native.evolve_model.domain());
        ARIADNE_TESTING_EXPECT_EQUAL(
            restored->reach_models.size(), native.reach_models.size());
        ARIADNE_TESTING_PROVE(step.representation_as<FlowStepModelType>() == nullptr);
    }

    Void test_reject_inconsistent_basic_result() {
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(
                StepSizeType(0u), UpperBoxType({{-1.0_x, 1.0_x}}),
                UpperBoxType({{-1.0_x, 1.0_x}})));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationStep(
                StepSizeType(0.03125_x), UpperBoxType({{-1.0_x, 1.0_x}}),
                UpperBoxType({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}})));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_validated_integration_step");
    ARIADNE_TESTING_ADD_TEST(TestValidatedIntegrationStep, test_single_valued_flow_preserves_patch);
    ARIADNE_TESTING_ADD_TEST(TestValidatedIntegrationStep, test_piecewise_reach_sections_cover_every_time);
    ARIADNE_TESTING_ADD_TEST(TestValidatedIntegrationStep, test_set_valued_box_step_without_forced_structure);
    ARIADNE_TESTING_ADD_TEST(TestValidatedIntegrationStep, test_structured_inclusion_retains_parameter_maps);
    ARIADNE_TESTING_ADD_TEST(TestValidatedIntegrationStep, test_reject_inconsistent_basic_result);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
