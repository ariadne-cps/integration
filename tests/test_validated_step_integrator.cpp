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
// Compile the shared step contract through the public module entry point.
#include "ariadne-integration.hpp"
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
        ARIADNE_TESTING_EXPECT_THROWS(std::logic_error,
            second.reach_sections(second.final_time()));
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
        // Autonomous flows keep a local [0,h] patch time although the
        // physical second step starts at h1.
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_offset(), h1);
        const auto local_sections = second.reach_sections(second.final_time());
        ARIADNE_TESTING_EXPECT_EQUAL(local_sections.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            local_sections[0u].domain(), second.evolve_model().domain());
        ARIADNE_TESTING_EXPECT_LESS(
            norm(local_sections[0u] - second.evolve_model()), 1e-8_pr);
        ARIADNE_TESTING_EXPECT_THROWS(std::out_of_range,
            second.reach_sections(StepSizeType(0u)));
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().domain(), initial);
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().argument_size(), initial.dimension());
        // The second tube begins at t=h1, not t=0. For x'=x with
        // x0 in [7/8,9/8] the exact second-step reach is
        // [7/8*exp(h1), 9/8*exp(h1+h2)] =
        // approximately [0.90277548, 1.17898988]. Test a strict
        // interior enclosure using dyadic endpoints (exactly representable
        // as ExactDouble), without assuming rounded exp().
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.90283203125_x, 1.178955078125_x}}), cast_exact_box(second.reach()))));
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
        // Use dyadic endpoints so the _x ExactDouble literals are exact.
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.15869140625_x, 0.15869140625_x}}), cast_exact_box(second.evolve()))));
    }

    Void test_nonautonomous_equation_continuation_uses_absolute_time() {
        // x'=t, x0 in [0,1/8], first step begins at t=1/2.
        // Exact final map: x0 + ((19/32)^2-(1/2)^2)/2
        //                = x0 + 105/2048.
        const auto t = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        const EffectiveVectorMultivariateFunction field = {t};
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.0625_x);
        const StepSizeType h2(0.03125_x);
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter adapter(taylor);
        const ValidatedStepIntegratorInterface& solver = adapter;
        const auto first = solver.step(ValidatedIntegrationRequest(
            field, initial, t0, h1, ExactBoxType(0u), bound));
        const auto second = solver.step(ValidatedIntegrationRequest(
            field, first, h2, bound));

        const auto x0 = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const EffectiveVectorMultivariateFunction expected = {
            x0 + Real(0.05126953125_x)};
        ARIADNE_TESTING_EXPECT_EQUAL(first.initial_time(), t0);
        ARIADNE_TESTING_EXPECT_EQUAL(first.final_time(), StepSizeType(0.5625_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), StepSizeType(0.5625_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.59375_x));
        ARIADNE_TESTING_PROVE(second.equation_has_explicit_time());
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_offset(), StepSizeType(0u));
        const auto absolute_sections = second.reach_sections(second.final_time());
        ARIADNE_TESTING_EXPECT_EQUAL(absolute_sections.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            absolute_sections[0u].domain(), second.evolve_model().domain());
        ARIADNE_TESTING_EXPECT_LESS(
            norm(absolute_sections[0u] - second.evolve_model()), 1e-8_pr);
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().domain(), initial);
        ARIADNE_TESTING_EXPECT_LESS(norm(second.evolve_model() - expected), 1e-4_pr);
    }

    Void test_constant_equation_parameters_are_not_resampled() {
        // x'=t+p, t0=1/2, p in [1/2,3/4]. p is selected once,
        // not independently at each step. The final map must keep
        // the ORIGINAL parameter coordinate and its entire domain.
        const auto t = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction field = {t + p};
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const ExactBoxType parameters({{0.5_x, 0.75_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.0625_x);
        const StepSizeType h2(0.03125_x);
        TaylorSeriesIntegrator taylor(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter adapter(taylor);
        const ValidatedStepIntegratorInterface& solver = adapter;

        const auto first = solver.step(ValidatedIntegrationRequest(
            field, initial, t0, h1, parameters, bound));
        const auto second = solver.step(ValidatedIntegrationRequest(
            field, first, h2, bound));

        const auto x0 = EffectiveScalarMultivariateFunction::coordinate(2u, 0u);
        const auto p0 = EffectiveScalarMultivariateFunction::coordinate(2u, 1u);
        const EffectiveVectorMultivariateFunction expected = {
            x0 + p0 * Real(0.09375_x) + Real(0.05126953125_x)};
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.equation_has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(second.constant_parameter_domain(), parameters);
        ARIADNE_TESTING_EXPECT_EQUAL(second.evolve_model().argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain(), product(initial, parameters));
        ARIADNE_TESTING_PROVE(second.reach_time_index().has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_index().value(), 2u);
        ARIADNE_TESTING_EXPECT_LESS(norm(second.evolve_model() - expected), 1e-4_pr);
    }

    Void test_nonautonomous_inclusion_fixed_parameters_and_switching_inputs() {
        // x'=t+p+u(t), p chosen once in [1/8,1/4] and
        // u(t) independently selectable from [-1/8,1/8].
        // From t=1/2 over h1=1/16,h2=1/32, exact final reachable set:
        // [105/2048,177/2048] = [0.05126953125,0.08642578125].
        // The box solver may be wider; it must enclose every trajectory.
        const auto t = EffectiveScalarMultivariateFunction::coordinate(4u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(4u, 2u);
        const auto u = EffectiveScalarMultivariateFunction::coordinate(4u, 3u);
        const EffectiveVectorMultivariateFunction field_function = {t + p + u};
        ValidatedTimedInputFieldBoxAdapterDP field(
            field_function, BoxDomainType({{-1.0_x, 1.0_x}}),
            IntervalDomainType(0.5_x, 0.75_x),
            BoxDomainType({{0.125_x, 0.25_x}}),
            BoxDomainType({{-0.125_x, 0.125_x}}));
        const ExactBoxType initial({{0.0_x, 0.0_x}});
        const ExactBoxType parameters({{0.125_x, 0.25_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.0625_x);
        const StepSizeType h2(0.03125_x);
        const FirstOrderMultifunctionIntegrator numerical_solver;
        const MultiflowStepIntegratorAdapter adapter(numerical_solver);
        const ValidatedStepIntegratorInterface& solver = adapter;

        const auto first = solver.step(ValidatedIntegrationRequest(
            field, initial, t0, h1, parameters, bound));
        const auto second = solver.step(ValidatedIntegrationRequest(
            field, first, h2, bound));

        ARIADNE_TESTING_EXPECT_EQUAL(first.initial_time(), t0);
        ARIADNE_TESTING_EXPECT_EQUAL(first.final_time(), StepSizeType(0.5625_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), StepSizeType(0.5625_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.59375_x));
        ARIADNE_TESTING_PROVE(second.has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(second.constant_parameter_domain(), parameters);
        ARIADNE_TESTING_PROVE(not second.has_evolve_model());
        ARIADNE_TESTING_PROVE(not second.has_reach_models());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.05126953125_x, 0.08642578125_x}}),
            cast_exact_box(second.evolve()))));
    }

    Void test_timed_inclusion_rejects_unsupported_capabilities_and_domains() {
        const auto t = EffectiveScalarMultivariateFunction::coordinate(4u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(4u, 2u);
        const auto u = EffectiveScalarMultivariateFunction::coordinate(4u, 3u);
        ValidatedTimedInputFieldBoxAdapterDP field(
            EffectiveVectorMultivariateFunction({t + p + u}),
            BoxDomainType({{-1.0_x, 1.0_x}}),
            IntervalDomainType(0.5_x, 0.75_x),
            BoxDomainType({{0.125_x, 0.25_x}}),
            BoxDomainType({{-0.125_x, 0.125_x}}));
        const ExactBoxType initial({{0.0_x, 0.0_x}});
        const ExactBoxType parameters({{0.125_x, 0.25_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        FirstOrderMultifunctionIntegrator first_order;
        MultiflowStepIntegratorAdapter first_order_adapter(first_order);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            first_order_adapter.step(ValidatedIntegrationRequest(
                field, initial, StepSizeType(0.71875_x),
                StepSizeType(0.0625_x), parameters, bound)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            first_order_adapter.step(ValidatedIntegrationRequest(
                field, initial, StepSizeType(0.5_x),
                StepSizeType(0.0625_x),
                ExactBoxType({{0.5_x, 0.75_x}}), bound)));
        ARIADNE_TESTING_EXPECT_THROWS(FlowBoundsException,
            first_order_adapter.step(ValidatedIntegrationRequest(
                field, initial, StepSizeType(0.5_x),
                StepSizeType(0.0625_x), parameters,
                UpperBoxType({{-0.0078125_x, 0.0078125_x}}))));

        // The high-order legacy solver has not been extended to time-dependent
        // inclusions, and must not be silently substituted by first order.
        TaylorSeriesIntegrator taylor(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                      lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high(
            taylor, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        MultiflowStepIntegratorAdapter unsupported(high);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            unsupported.step(ValidatedIntegrationRequest(
                field, initial, StepSizeType(0.5_x),
                StepSizeType(0.0625_x), parameters, bound)));
    }

    Void test_timed_high_order_inclusion_preserves_constant_parameters_and_inputs() {
        // Lift (x,t,p)'=(t+p+u(t),1,0) to the historical autonomous
        // inclusion solver. After two steps, old input coefficients
        // and the fixed parameter p must retain their original domains.
        const auto t = EffectiveScalarMultivariateFunction::coordinate(4u, 1u);
        const auto p = EffectiveScalarMultivariateFunction::coordinate(4u, 2u);
        const auto u = EffectiveScalarMultivariateFunction::coordinate(4u, 3u);
        const EffectiveVectorMultivariateFunction f = {t + p + u};
        ValidatedTimedInputFieldBoxAdapterDP field(
            f, BoxDomainType({{-1.0_x, 1.0_x}}),
            IntervalDomainType(0.0_x, 1.0_x),
            BoxDomainType({{0.0_x, 1.0_x}}),
            BoxDomainType({{-0.125_x, 0.125_x}}));
        const ExactBoxType initial({{0.0_x, 0.0_x}});
        const ExactBoxType parameters({{0.125_x, 0.25_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);

        TaylorSeriesIntegrator taylor(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                      lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high(
            taylor, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter adapter(high);
        const ValidatedStepIntegratorInterface& integrator = adapter;
        const auto first = integrator.step(ValidatedIntegrationRequest(
            field, initial, t0, h1, parameters, bound));
        const auto second = integrator.step(ValidatedIntegrationRequest(
            field, first, h2, bound));

        ARIADNE_TESTING_PROVE(first.has_evolve_model());
        ARIADNE_TESTING_PROVE(first.has_reach_models());
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), StepSizeType(0.53125_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.546875_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.constant_parameter_domain(), parameters);
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain()[range(0u, 2u)], product(initial, parameters));
        ARIADNE_TESTING_PROVE(
            second.evolve_model().argument_size() > first.evolve_model().argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain()[range(0u, first.evolve_model().argument_size())],
            first.evolve_model().domain());
        ARIADNE_TESTING_PROVE(second.reach_time_index().has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_index().value(), 2u);
        for (const auto& reach_patch: second.reach_models()) {
            ARIADNE_TESTING_EXPECT_EQUAL(
                reach_patch.domain()[2u],
                IntervalDomainType(0.53125_x, 0.546875_x));
        }
        // The lifted inclusion solver already stores absolute time, even
        // when the time coordinate follows original fixed parameters.
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_offset(), StepSizeType(0u));
        const auto inclusion_sections = second.reach_sections(second.final_time());
        ARIADNE_TESTING_PROVE(not inclusion_sections.empty());
        for (const auto& section: inclusion_sections) {
            // Reach uses ZeroApproximation and evolve ConstantApproximation.
            // Their newly introduced input coordinates need not agree:
            // reach may carry no fresh input coordinate, while evolve
            // introduces one. Both must retain the ORIGINAL initial-state,
            // fixed-parameter and earlier input-uncertainty coordinates.
            const SizeType retained = first.evolve_model().argument_size();
            ARIADNE_TESTING_PROVE(section.argument_size() >= retained);
            ARIADNE_TESTING_EXPECT_EQUAL(
                section.domain()[range(0u, retained)],
                first.evolve_model().domain());
            ARIADNE_TESTING_EXPECT_EQUAL(
                section.domain()[range(0u, retained)],
                second.evolve_model().domain()[range(0u, retained)]);
        }

        // A fixed parameter is selected ONCE for the whole trajectory:
        // after H=1/32+1/64=3/64, changing p from 1/8 to 1/4
        // must change the endpoint by H*(1/4-1/8)=3/512.
        // Slice the SAME validated endpoint map while retaining all
        // original state and independently selected input coordinates.
        const auto fixed_low = partial_evaluate(second.evolve_model(), 1u, 0.125_x);
        const auto fixed_high = partial_evaluate(second.evolve_model(), 1u, 0.25_x);
        ARIADNE_TESTING_EXPECT_EQUAL(fixed_low.domain(), fixed_high.domain());
        const auto parameter_effect = fixed_high - fixed_low;
        ARIADNE_TESTING_NOTE("fixed-parameter endpoint difference="
            << parameter_effect.range() << "; expected=" << 0.005859375_x);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.005859375_x, 0.005859375_x}}),
            cast_exact_box(UpperBoxType(parameter_effect.range())))));

        // Exact first endpoint: [33/2048,57/2048].
        // Exact second endpoint: [201/8192,345/8192].
        // Keep the full analytical enclosure assertion below; these
        // recorded values localize failures without weakening coverage.
        ARIADNE_TESTING_NOTE("timed high-order first endpoint box=" << first.evolve()
            << "; first model range=" << first.evolve_model().range());
        ARIADNE_TESTING_NOTE("timed high-order second endpoint box=" << second.evolve()
            << "; second model range=" << second.evolve_model().range());
        ARIADNE_TESTING_NOTE("time t0=" << t0 << "; first h=" << h1
            << "; second h=" << h2 << "; parameters=" << parameters);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.01611328125_x, 0.02783203125_x}}),
            cast_exact_box(first.evolve()))));
        // Check the two exact extrema separately to distinguish missing
        // low trajectories from missing high trajectories.
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0245361328125_x, 0.0245361328125_x}}),
            cast_exact_box(second.evolve()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0421142578125_x, 0.0421142578125_x}}),
            cast_exact_box(second.evolve()))));
        // Exact endpoint x'=t+p+u(t), |u|<=1/8, p in[1/8,1/4]
        // after total H=3/64: [201/8192,345/8192].
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0245361328125_x, 0.0421142578125_x}}),
            cast_exact_box(second.evolve()))));
    }

    Void test_timed_high_order_rejects_missing_structural_capability() {
        // Box evaluation alone cannot justify a guessed high-order
        // representation or selecting a different integrator.
        class OnlyTimedBoxes final : public ValidatedTimedBoxMultifunctionInterface {
          public:
            SizeType argument_size() const override { return 1u; }
            SizeType result_size() const override { return 1u; }
            BoxDomainType domain() const override {
                return BoxDomainType({{-1.0_x, 1.0_x}});
            }
            IntervalDomainType time_domain() const override {
                return IntervalDomainType(0.0_x, 1.0_x);
            }
            BoxDomainType parameter_domain() const override {
                return BoxDomainType({{0.0_x, 1.0_x}});
            }
            UpperBoxType evaluate(const UpperBoxType&, const UpperIntervalType&,
                                  const UpperBoxType&) const override {
                return UpperBoxType({{-1.0_x, 1.0_x}});
            }
        } field;
        TaylorSeriesIntegrator taylor(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                      lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high(
            taylor, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter adapter(high);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            adapter.step(ValidatedIntegrationRequest(
                field, ExactBoxType({{0.0_x, 0.0_x}}), StepSizeType(0.5_x),
                StepSizeType(0.03125_x), ExactBoxType({{0.125_x, 0.25_x}}),
                UpperBoxType({{-1.0_x, 1.0_x}}))));
    }

    Void test_equation_then_timed_inclusion_preserves_constant_parameter() {
        // First use a native nonautonomous equation solver for
        // x'=t+p. Then use the autonomous high-order inclusion solver
        // (internally lifted) for x'=t+p+u(t). The original p coordinate
        // must be shared across both solvers; only u(t) is a fresh input.
        const auto t_eq = EffectiveScalarMultivariateFunction::coordinate(3u, 1u);
        const auto p_eq = EffectiveScalarMultivariateFunction::coordinate(3u, 2u);
        const EffectiveVectorMultivariateFunction equation = {t_eq + p_eq};

        const auto t_inc = EffectiveScalarMultivariateFunction::coordinate(4u, 1u);
        const auto p_inc = EffectiveScalarMultivariateFunction::coordinate(4u, 2u);
        const auto u_inc = EffectiveScalarMultivariateFunction::coordinate(4u, 3u);
        const EffectiveVectorMultivariateFunction inclusion = {t_inc + p_inc + u_inc};
        ValidatedTimedInputFieldBoxAdapterDP set_valued(
            inclusion, BoxDomainType({{-1.0_x, 1.0_x}}),
            IntervalDomainType(0.0_x, 1.0_x),
            BoxDomainType({{0.0_x, 1.0_x}}),
            BoxDomainType({{-0.125_x, 0.125_x}}));

        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const ExactBoxType parameters({{0.125_x, 0.25_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType t0(0.5_x);
        const StepSizeType h1(0.03125_x);
        const StepSizeType h2(0.015625_x);

        TaylorSeriesIntegrator equation_solver(step_maximum_error = 1e-6, order = 8);
        const FlowStepIntegratorAdapter equation_adapter(equation_solver);
        const ValidatedStepIntegratorInterface& first_solver = equation_adapter;
        const auto first = first_solver.step(ValidatedIntegrationRequest(
            equation, initial, t0, h1, parameters, bound));

        TaylorSeriesIntegrator high_taylor(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                           lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high_inclusion(
            high_taylor, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const StructuredInclusionStepIntegratorAdapter inclusion_adapter(high_inclusion);
        const ValidatedStepIntegratorInterface& second_solver = inclusion_adapter;
        const auto second = second_solver.step(ValidatedIntegrationRequest(
            set_valued, first, h2, bound));

        ARIADNE_TESTING_PROVE(first.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.has_evolve_model());
        ARIADNE_TESTING_PROVE(second.has_reach_models());
        ARIADNE_TESTING_EXPECT_EQUAL(second.initial_time(), StepSizeType(0.53125_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.final_time(), StepSizeType(0.546875_x));
        ARIADNE_TESTING_EXPECT_EQUAL(second.constant_parameter_domain(), parameters);
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain()[range(0u, 2u)], product(initial, parameters));
        ARIADNE_TESTING_PROVE(
            second.evolve_model().argument_size() > first.evolve_model().argument_size());
        ARIADNE_TESTING_EXPECT_EQUAL(
            second.evolve_model().domain()[range(0u, first.evolve_model().argument_size())],
            first.evolve_model().domain());
        ARIADNE_TESTING_PROVE(second.reach_time_index().has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(second.reach_time_index().value(), 2u);

        // The exact endpoint has H=3/64; u is active only during h2=1/64.
        // x0 in [0,1/8], p in [1/8,1/4], u(t) in [-1/8,1/8]:
        // [H*t0+H^2/2+H/8-h2/8,
        //  H*t0+H^2/2+H/4+h2/8+1/8]
        // = [233/8192, 1337/8192].
        ARIADNE_TESTING_NOTE("equation-to-inclusion first=" << first.evolve()
            << "; second=" << second.evolve());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0284423828125_x, 0.1632080078125_x}}),
            cast_exact_box(second.evolve()))));

        // Compare two slices of the SAME endpoint map. The difference
        // caused by p must be H*(1/4-1/8)=3/512 regardless of the
        // algorithm used in the first step.
        const auto low_p = partial_evaluate(second.evolve_model(), 1u, 0.125_x);
        const auto high_p = partial_evaluate(second.evolve_model(), 1u, 0.25_x);
        ARIADNE_TESTING_EXPECT_EQUAL(low_p.domain(), high_p.domain());
        const auto effect = high_p - low_p;
        ARIADNE_TESTING_NOTE("equation-to-inclusion parameter effect=" << effect.range());
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.005859375_x, 0.005859375_x}}),
            cast_exact_box(UpperBoxType(effect.range())))));
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
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_nonautonomous_equation_continuation_uses_absolute_time);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_constant_equation_parameters_are_not_resampled);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_nonautonomous_inclusion_fixed_parameters_and_switching_inputs);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_timed_inclusion_rejects_unsupported_capabilities_and_domains);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_timed_high_order_inclusion_preserves_constant_parameters_and_inputs);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_timed_high_order_rejects_missing_structural_capability);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_equation_then_timed_inclusion_preserves_constant_parameter);
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_reject_invalid_request);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
