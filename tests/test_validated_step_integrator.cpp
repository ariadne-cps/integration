/***************************************************************************
 *            test_validated_step_integrator.cpp
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

/** @file test_validated_step_integrator.cpp
 *  @brief Exercise one polymorphic integration interface with distinct fields.
 */

#include <stdexcept>

#include "function/symbolic_function.hpp"
#include "integration/validated_step_integrator.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestValidatedStepIntegrator : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_common_polymorphic_step_for_both_field_capabilities() {
        const ValidatedVectorMultivariateFunction ordinary = {
            ValidatedScalarMultivariateFunction::constant(1u, 1)};
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction input_field =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP set_valued(input_field,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));

        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        const ValidatedIntegrationRequest function_request(ordinary, initial, h, bound);
        const ValidatedIntegrationRequest multifunction_request(set_valued, initial, h, bound);

        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        FirstOrderMultifunctionIntegrator low_order;
        const FlowStepIntegratorAdapter flow_solver(taylor);
        const MultiflowStepIntegratorAdapter inclusion_solver(low_order);

        // The caller only depends on this one numerical step interface.
        const ValidatedStepIntegratorInterface* const solvers[2] = {&flow_solver, &inclusion_solver};
        const auto ordinary_step = solvers[0]->step(function_request);
        const auto inclusion_step = solvers[1]->step(multifunction_request);

        ARIADNE_TESTING_EXPECT_EQUAL(ordinary_step.time_step(), h);
        ARIADNE_TESTING_EXPECT_EQUAL(inclusion_step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.03125_x, 0.15625_x}}),
            cast_exact_box(ordinary_step.evolve()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}),
            cast_exact_box(inclusion_step.evolve()))));
        ARIADNE_TESTING_PROVE(ordinary_step.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(
            ordinary_step.evolve_model().argument_size(), initial.dimension());
        ARIADNE_TESTING_PROVE(ordinary_step.has_reach_models());
        ARIADNE_TESTING_PROVE(not inclusion_step.has_evolve_model());
        ARIADNE_TESTING_PROVE(not inclusion_step.has_reach_models());
        ARIADNE_TESTING_PROVE(
            ordinary_step.representation_as<FlowStepModelType>() != nullptr);
        ARIADNE_TESTING_PROVE(not inclusion_step.has_representation());

        // An unsuitable numerical input must not select another solver.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            solvers[0]->step(multifunction_request));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            solvers[1]->step(function_request));
    }

    Void test_structured_inclusion_retains_native_maps() {
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP field(f,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator selected(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter solver(selected);
        const ValidatedStepIntegratorInterface& common = solver;
        const auto step = common.step(ValidatedIntegrationRequest(field, initial, h, bound));

        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}),
            cast_exact_box(step.evolve()))));
        // The future evolver need not inspect the concrete solver result.
        ARIADNE_TESTING_PROVE(step.has_reach_models());
        ARIADNE_TESTING_PROVE(step.has_evolve_model());
        ARIADNE_TESTING_PROVE(step.reach_models().size() >= 1u);
        ARIADNE_TESTING_PROVE(step.evolve_model().argument_size() > initial.dimension());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}),
            cast_exact_box(UpperBoxType(step.evolve_model().range())))));
        // The original algorithm-specific result remains available, if needed.
        const auto* original = step.representation_as<SelectedInclusionStep>();
        ARIADNE_TESTING_PROVE(original != nullptr);
        ARIADNE_TESTING_EXPECT_EQUAL(original->time_step, h);
    }

    Void test_structured_two_step_continuation_uses_common_result() {
        // x'=u(t), |u(t)|<=1. Each step must admit a fresh independent
        // input while preserving the coordinates and inputs of earlier steps.
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP field(f,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));
        const ExactBoxType initial({{-0.125_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator selected(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter adapter(selected);
        const ValidatedStepIntegratorInterface& solver = adapter;
        const auto first = solver.step(
            ValidatedIntegrationRequest(field, initial, h1, bound));
        const auto second = solver.step(
            ValidatedIntegrationRequest(field, first, h2, bound));

        ARIADNE_TESTING_EXPECT_EQUAL(first.initial_time(), StepSizeType(0u));
        ARIADNE_TESTING_EXPECT_EQUAL(first.final_time(), h1);
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), h1);
        ARIADNE_TESTING_EXPECT_EQUAL(second.time_step(), h2);
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.046875_x));
        ARIADNE_TESTING_PROVE(first.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_PROVE(
            second.evolve_model().argument_size() > first.evolve_model().argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain()[range(0u, first.evolve_model().argument_size())],
            first.evolve_model().domain());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.171875_x, 0.171875_x}}),
            cast_exact_box(second.evolve()))));
        ARIADNE_TESTING_PROVE(second.has_reach_models());
        for (const auto& reach_model: second.reach_models()) {
            ARIADNE_TESTING_EXPECT_EQUAL(
                reach_model.domain()[initial.dimension()],
                IntervalDomainType(0.03125_x, 0.046875_x));
        }
    }

    Void test_box_only_continuation_keeps_hull_and_time() {
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP field(f,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));
        const ExactBoxType initial({{-0.125_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);
        FirstOrderMultifunctionIntegrator low;
        const MultiflowStepIntegratorAdapter adapter(low);
        const ValidatedStepIntegratorInterface& solver = adapter;
        const auto first = solver.step(ValidatedIntegrationRequest(field, initial, h1, bound));
        const auto second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound));

        ARIADNE_TESTING_PROVE(not first.has_evolve_model());
        ARIADNE_TESTING_PROVE(not second.has_evolve_model());
        ARIADNE_TESTING_PROVE(not second.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), h1);
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.046875_x));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.171875_x, 0.171875_x}}),
            cast_exact_box(second.evolve()))));
    }

    Void test_no_implicit_loss_of_structured_continuation() {
        // A structured map must not be silently restarted as an endpoint box
        // by a solver which cannot preserve its accumulated input parameters.
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP field(f,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));
        const ExactBoxType initial({{-0.125_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator selected(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter high(selected);
        const auto first = high.step(ValidatedIntegrationRequest(field, initial, h, bound));
        FirstOrderMultifunctionIntegrator low;
        const MultiflowStepIntegratorAdapter box_adapter(low);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            box_adapter.step(ValidatedIntegrationRequest(field, first, h, bound)));

        // The concrete type of the preceding solver is not consulted.
        ARIADNE_TESTING_PROVE(first.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(
            ValidatedIntegrationRequest::continuation_domain(first), initial);
    }

    Void test_equation_continuation_preserves_validated_mapping() {
        // x'=x, two unequal steps: continuation must retain dependence on
        // the original initial coordinate, not replace the endpoint map by a box.
        const auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const EffectiveVectorMultivariateFunction field = {x};
        const ExactBoxType initial({{0.875_x, 1.125_x}});
        const UpperBoxType bound({{0.5_x, 2.0_x}});
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter adapter(taylor);
        const ValidatedStepIntegratorInterface& solver = adapter;

        const auto first = solver.step(ValidatedIntegrationRequest(field, initial, h1, bound));
        const auto second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound));
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), h1);
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.046875_x));
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_models().size(), 1u);
        ARIADNE_TESTING_PROVE(second.reach_time_index().has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_index().value(), initial.dimension());
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().domain(), initial);
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().argument_size(), initial.dimension());
        // The second tube begins at t=h1, not t=0. For x'=x with
        // x0 in [7/8,9/8] the exact second-step reach is
        // [7/8*exp(h1), 9/8*exp(h1+h2)] =
        // approximately [0.90277548, 1.17898988]. Test a strict
        // interior enclosure to avoid assuming exactly rounded exp().
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.9028_x, 1.1789_x}}), cast_exact_box(second.reach()))));
        // Analytical endpoint x0*exp(h1+h2), with a tolerance accounting
        // for the validated Taylor-model remainder and outward rounding.
        const EffectiveVectorMultivariateFunction expected = {
            x * exp(Real(0.046875_x))};
        ARIADNE_TESTING_EXPECT_LESS(norm(second.evolve_model() - expected), 1e-4_pr);
    }

    Void test_equation_continuation_retains_prior_input_parameters() {
        // An equation solver receives a structured endpoint map from a
        // preceding inclusion. Every prior input parameter must survive.
        const auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction varying_field =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        ValidatedInputFieldBoxAdapterDP inclusion_field(varying_field,
            BoxDomainType({{-1.0_x, 1.0_x}}),
            BoxDomainType({{-1.0_x, 1.0_x}}));
        const auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const EffectiveVectorMultivariateFunction equation_field = {x};
        const ExactBoxType initial({{-0.125_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);

        TaylorSeriesIntegrator taylor(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                      lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator inclusion(
            taylor, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter first_solver(inclusion);
        const FlowStepIntegratorAdapter second_solver(taylor);
        const auto first = first_solver.step(
            ValidatedIntegrationRequest(inclusion_field, initial, h1, bound));
        const auto second = second_solver.step(
            ValidatedIntegrationRequest(equation_field, first, h2, bound));

        ARIADNE_TESTING_PROVE(first.has_evolve_model());
        ARIADNE_TESTING_PROVE(first.evolve_model().argument_size() > initial.dimension());
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().argument_size(), first.evolve_model().argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().domain(), first.evolve_model().domain());
        ARIADNE_TESTING_PROVE(second.reach_time_index().has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.reach_time_index().value(), first.evolve_model().argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.046875_x));
        // The exact endpoint is [-5/32,5/32]*exp(1/64) =
        // approximately [-0.1587105795,0.1587105795].
        // The old [-0.1612,0.1612] expectation was wider than
        // the reachable set, and hence was not a valid lower bound.
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.1587_x, 0.1587_x}}), cast_exact_box(second.evolve()))));
    }

    Void test_reject_invalid_request() {
        const ValidatedVectorMultivariateFunction ordinary = {
            ValidatedScalarMultivariateFunction::constant(1u, 1)};
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationRequest(
                ordinary, ExactBoxType({{0.0_x, 0.0_x}}),
                StepSizeType(0u), bound));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedIntegrationRequest(
                ordinary, ExactBoxType({{0.0_x, 0.0_x}}),
                StepSizeType(0.03125_x),
                UpperBoxType({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}})));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_validated_step_integrator");
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_common_polymorphic_step_for_both_field_capabilities);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_structured_inclusion_retains_native_maps);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_structured_two_step_continuation_uses_common_result);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_box_only_continuation_keeps_hull_and_time);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_no_implicit_loss_of_structured_continuation);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_equation_continuation_preserves_validated_mapping);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_equation_continuation_retains_prior_input_parameters);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_reject_invalid_request);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
