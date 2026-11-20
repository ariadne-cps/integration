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


    Void test_every_endpoint_strategy_encloses_full_input_range() {
        // x' = u(t), u(t) in [-1,1], x(0) in [-1/8,1/8].
        // For any measurable input, the entire endpoint set at h=1/32 is
        // [-5/32,5/32]. No evolve strategy may freeze or discard inputs.
        auto u = EffectiveFormula::coordinate(1u);
        const EffectiveVectorMultivariateFunction f =
            EffectiveVectorMultivariateFormulaFunction(
                2u, Vector<EffectiveFormula>{u});
        const BoxDomainType U({{-1.0_x, 1.0_x}});
        const BoxDomainType initial_domain({{-0.125_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        const ExactBoxType analytical_endpoint({{-0.15625_x, 0.15625_x}});

        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        InclusionIntegratorFactory factory(base);
        const auto initial = base.function_factory().create_identity(initial_domain);
        const InputApproximation rigorous_reach(ZeroApproximation());
        for (const auto& endpoint: {InputApproximation(ZeroApproximation()),
                                   InputApproximation(ConstantApproximation()),
                                   InputApproximation(AffineApproximation()),
                                   InputApproximation(SinusoidalApproximation()),
                                   InputApproximation(PiecewiseApproximation())}) {
            SelectedInclusionIntegrator selected(
                factory, f, U, rigorous_reach, endpoint);
            ARIADNE_TESTING_EXPECT_EQUAL(selected.reach_index(), rigorous_reach.reference().index());
            ARIADNE_TESTING_EXPECT_EQUAL(selected.evolve_index(), endpoint.reference().index());
            const auto step = selected.step(
                initial_domain, initial, bound, TimeStepType(0u), h);
            const auto boxes = step.bounding_step(h);
            ARIADNE_TESTING_PROVE(definitely(subset(
                analytical_endpoint, cast_exact_box(boxes.evolve))));
            ARIADNE_TESTING_PROVE(definitely(subset(
                analytical_endpoint, cast_exact_box(boxes.reach))));
        }
    }

    Void test_parametric_box_extension_and_low_order_bridge() {
        BoxDomainType state_domain({{-1.0_x, 1.0_x}});
        // The owned structural representation outlives these temporaries.
        ValidatedInputFieldBoxAdapterDP field_view(field(), inputs(), state_domain);
        UpperBoxType evaluation_box({{0.0_x, 0.125_x}});
        auto range = field_view.evaluate(evaluation_box);
        ARIADNE_TESTING_EXPECT_EQUAL(range.dimension(), 1u);
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{-0.125_x, 0.25_x}}), cast_exact_box(range))));

        FirstOrderMultifunctionIntegrator low_order;
        const MultifunctionIntegratorInterface& common = low_order;
        auto box_step = common.flow_step(field_view, ExactBoxType({{0.0_x, 0.125_x}}),
                                         StepSizeType(0.03125_x), UpperBoxType({{-1.0_x, 1.0_x}}));
        ARIADNE_TESTING_EXPECT_EQUAL(box_step.time_step, StepSizeType(0.03125_x));
        ARIADNE_TESTING_PROVE(definitely(inside(box_step.reach, ExactBoxType({{-1.0_x, 1.0_x}}))));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(box_step.evolve)));

        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            field_view.evaluate(UpperBoxType({{-2.0_x, 2.0_x}})));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedInputFieldBoxAdapterDP(field(), BoxDomainType(0u), state_domain));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedInputFieldBoxAdapterDP(field(), inputs(), BoxDomainType({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}})));
    }

    Void test_high_order_common_integrator_bridge() {
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high_order(base,
            InputApproximation(ZeroApproximation()),
            InputApproximation(AffineApproximation()));
        ARIADNE_TESTING_EXPECT_EQUAL(high_order.reach_index(), ZeroApproximation().index());
        ARIADNE_TESTING_EXPECT_EQUAL(high_order.evolve_index(), AffineApproximation().index());

        ValidatedInputFieldBoxAdapterDP structural(field(), inputs(), BoxDomainType({{-1.0_x, 1.0_x}}));
        ExactBoxType initial({{0.0_x, 0.125_x}});
        UpperBoxType bound({{-1.0_x, 1.0_x}});
        StepSizeType h(0.03125_x);
        const MultifunctionIntegratorInterface& common = high_order;
        const auto box_step = common.flow_step(structural, initial, h, bound);
        ARIADNE_TESTING_EXPECT_EQUAL(box_step.time_step, h);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(box_step.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(box_step.evolve)));
        const auto patch_step = high_order.structured_step(structural, initial, h, bound);
        ARIADNE_TESTING_EXPECT_EQUAL(patch_step.reach_models.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(patch_step.evolve_model.argument_size(), 3u);
        ARIADNE_TESTING_PROVE(same(cast_exact_box(patch_step.bounding_step(h).reach),
                                  cast_exact_box(box_step.reach)));

        const auto expected = FloatDPBounds(0.0625_x, dp) * exp(FloatDPBounds(h, dp));
        UpperIntervalType nominal(expected.lower(), expected.upper());
        ARIADNE_TESTING_PROVE(refines(nominal, box_step.evolve[0u]));

        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            high_order.structured_step(structural, initial, h, UpperBoxType({{-2.0_x, 2.0_x}})));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            SelectedParametricMultifunctionIntegrator(base,
                InputApproximation(AffineApproximation()),
                InputApproximation(ConstantApproximation())));
    }

    Void test_high_order_interval_taylor_via_parametric_conversion() {
        TaylorSeriesIntegrator base(step_maximum_error = 1e-6, order = 8);
        SelectedParametricMultifunctionIntegrator high_order(base,
            InputApproximation(ConstantApproximation()), InputApproximation(ConstantApproximation()));
        BoxDomainType domain({{-1.0_x, 1.0_x}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        using Model = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
        Scalar uncertain = Scalar::constant(domain, Interval<FloatDPUpperBound>(-1, 1, dp), sweeper);
        ValidatedIntervalTaylorBoxMultifunctionDP unstructured(Model({uncertain}));
        const MultifunctionIntegratorInterface& common = high_order;
        const auto representation = unstructured.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(representation.field.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(representation.field.argument_size(),
                                     1u + representation.inputs.dimension());
        ARIADNE_TESTING_PROVE(representation.inputs.dimension() > 0u);

        // The same high-order integrator now accepts the interval Taylor
        // model through the capability interface rather than concrete RTTI.
        auto step = common.flow_step(unstructured, ExactBoxType({{0.0_x, 0.125_x}}),
                                     StepSizeType(0.03125_x), UpperBoxType({{-1.0_x, 1.0_x}}));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolve)));
    }

    Void test_interval_taylor_conversion_of_state_dependent_input() {
        // x' = x*u(t), u(t) in [-1,1] as an interval-coefficient Taylor model.
        BoxDomainType domain({{0.5_x, 1.5_x}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        using Model = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
        Scalar x = Scalar::coordinate(domain, 0u, sweeper);
        Scalar field = Interval<FloatDPUpperBound>(-1, 1, dp) * x;
        ValidatedIntervalTaylorBoxMultifunctionDP vector_field(Model({field}));

        // Check model evaluation before any conversion is requested.
        Vector<FloatDPBounds> scalar_argument({FloatDPBounds(1.0_x, dp)});
        const auto scalar_image = field(scalar_argument);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(scalar_image)));
        const UpperBoxType early_states({{0.875_x, 1.125_x}});
        const auto early_image = vector_field.evaluate(early_states);
        ARIADNE_TESTING_EXPECT_EQUAL(early_image.dimension(), 1u);

        const auto form = vector_field.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(form.field.result_size(), 1u);
        ARIADNE_TESTING_PROVE(form.inputs.dimension() > 0u);

        // For all x and all uncertain coefficients, the converted field
        // must contain the original interval Taylor field value.
        const UpperBoxType states({{0.875_x, 1.125_x}});
        const auto interval_image = vector_field.evaluate(states);
        ARIADNE_TESTING_EXPECT_EQUAL(interval_image.dimension(), 1u);
        const Vector<FloatDPBounds> arguments(
            1u + form.inputs.dimension(), [&](SizeType j) {
                if (j == 0u) {
                    return FloatDPBounds(states[0u].lower_bound(), states[0u].upper_bound());
                }
                const auto& input = form.inputs[j - 1u];
                return FloatDPBounds(input.lower_bound(), input.upper_bound());
            });
        const Vector<FloatDPBounds> expanded_image = form.field(arguments);
        ARIADNE_TESTING_EXPECT_EQUAL(expanded_image.size(), 1u);
        const UpperBoxType expanded(
            expanded_image.size(), [&](SizeType j) {
                return UpperIntervalType(expanded_image[j].lower(), expanded_image[j].upper());
            });
        ARIADNE_TESTING_PROVE(definitely(subset(interval_image, cast_exact_box(expanded))));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(expanded)));

        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator integrator(base,
            InputApproximation(ZeroApproximation()), InputApproximation(ConstantApproximation()));
        ARIADNE_TESTING_EXPECT_EQUAL(integrator.reach_index(), ZeroApproximation().index());
        const auto step = integrator.flow_step(vector_field, ExactBoxType({{1.0_x, 1.0_x}}),
                                               StepSizeType(0.03125_x), UpperBoxType({{0.5_x, 1.5_x}}));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolve)));
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{1.0_x, 1.0_x}}),
                                                 cast_exact_box(step.reach))));
        // x(t) = exp(integral_0^t u(s) ds), |u(s)| <= 1:
        // all endpoints in [exp(-h), exp(h)] must be enclosed, not just x(0).
        const FloatDPBounds h_bounds(StepSizeType(0.03125_x), dp);
        const FloatDPBounds minimum = exp(-h_bounds);
        const FloatDPBounds maximum = exp(h_bounds);
        const UpperIntervalType exact_endpoint_hull(minimum.lower(), maximum.upper());
        ARIADNE_TESTING_PROVE(refines(exact_endpoint_hull, step.evolve[0u]));
        ARIADNE_TESTING_PROVE(refines(exact_endpoint_hull, step.reach[0u]));
    }

    Void test_interval_taylor_conversion_preserves_switching_endpoints() {
        // x' = u(t), y' = x with u(t) in [-1,1].
        // Switching u=+1 on [0,h/2] and u=-1 on [h/2,h]
        // reaches x(h)=0, y(h)=h*h/4. A frozen coefficient cannot do this.
        BoxDomainType domain({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        using Model = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
        Scalar uncertain = Scalar::constant(domain, Interval<FloatDPUpperBound>(-1, 1, dp), sweeper);
        Scalar x = Scalar::coordinate(domain, 0u, sweeper);
        ValidatedIntervalTaylorBoxMultifunctionDP field(Model({uncertain, x}));

        const auto form = field.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(form.field.result_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(form.field.argument_size(), 2u + form.inputs.dimension());

        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator integrator(base,
            InputApproximation(ZeroApproximation()), InputApproximation(AffineApproximation()));
        const StepSizeType h(0.03125_x);
        const auto step = integrator.flow_step(field,
            ExactBoxType({{0.0_x, 0.0_x}, {0.0_x, 0.0_x}}), h,
            UpperBoxType({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}}));

        const ExactBoxType positive_constant_input({
            {0.03125_x, 0.03125_x},
            {0.00048828125_x, 0.00048828125_x}});
        const ExactBoxType negative_constant_input({
            {-0.03125_x, -0.03125_x},
            {-0.00048828125_x, -0.00048828125_x}});
        const ExactBoxType switching_input({
            {0.0_x, 0.0_x},
            {0.000244140625_x, 0.000244140625_x}});

        ARIADNE_TESTING_PROVE(definitely(subset(positive_constant_input, cast_exact_box(step.evolve))));
        ARIADNE_TESTING_PROVE(definitely(subset(negative_constant_input, cast_exact_box(step.evolve))));
        ARIADNE_TESTING_PROVE(definitely(subset(switching_input, cast_exact_box(step.evolve))));
        ARIADNE_TESTING_PROVE(definitely(subset(switching_input, cast_exact_box(step.reach))));
    }


    Void test_multidimensional_interval_taylor_conversion() {
        // The converted parametric field must enclose the interval Taylor
        // model on boxes with both signs and mixed state dependencies.
        BoxDomainType domain({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        using Model = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
        Scalar x = Scalar::coordinate(domain, 0u, sweeper);
        Scalar y = Scalar::coordinate(domain, 1u, sweeper);
        Scalar first = Interval<FloatDPUpperBound>(-1, 1, dp) * x +
                       Interval<FloatDPUpperBound>(0, 2, dp) * y;
        Scalar second = Interval<FloatDPUpperBound>(-1, 1, dp) * y - x;
        ValidatedIntervalTaylorBoxMultifunctionDP field(Model({first, second}));
        const auto form = field.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(form.field.result_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(form.field.argument_size(),
                                     2u + form.inputs.dimension());

        for (const auto& states: {
                UpperBoxType({{0.5_x, 0.75_x}, {0.25_x, 0.5_x}}),
                UpperBoxType({{-0.75_x, -0.5_x}, {-0.5_x, -0.25_x}}),
                UpperBoxType({{-0.5_x, 0.5_x}, {-0.5_x, 0.5_x}}),
                UpperBoxType({{0.0_x, 0.0_x}, {0.0_x, 0.0_x}})}) {
            const auto original = field.evaluate(states);
            const Vector<FloatDPBounds> arguments(
                2u + form.inputs.dimension(), [&](SizeType j) {
                    if (j < 2u) {
                        return FloatDPBounds(states[j].lower_bound(), states[j].upper_bound());
                    }
                    const auto& input = form.inputs[j - 2u];
                    return FloatDPBounds(input.lower_bound(), input.upper_bound());
                });
            const Vector<FloatDPBounds> expanded_values = form.field(arguments);
            const UpperBoxType expanded(expanded_values.size(), [&](SizeType j) {
                return UpperIntervalType(expanded_values[j].lower(), expanded_values[j].upper());
            });
            ARIADNE_TESTING_EXPECT_EQUAL(original.dimension(), 2u);
            ARIADNE_TESTING_EXPECT_EQUAL(expanded.dimension(), 2u);
            ARIADNE_TESTING_PROVE(definitely(subset(original, cast_exact_box(expanded))));
        }
    }

    Void test_shared_input_correlation_is_explicit() {
        // The structured field (u,-u) has ONE shared time-varying input.
        // A pair of interval Taylor coefficient ranges alone cannot encode
        // this anticorrelation; conversion must conservatively keep the
        // component coefficients independent instead of inventing a link.
        const BoxDomainType domain({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}});
        const BoxDomainType input_domain({{-1.0_x, 1.0_x}});
        const auto u = EffectiveFormula::coordinate(2u);
        const EffectiveVectorMultivariateFunction correlated =
            EffectiveVectorMultivariateFormulaFunction(
                3u, Vector<EffectiveFormula>{u, -u});
        ValidatedInputFieldBoxAdapterDP structured(correlated, input_domain, domain);
        const auto shared = structured.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(shared.inputs.dimension(), 1u);
        const Vector<FloatDPBounds> same_input({
            FloatDPBounds(0.0_x, dp),
            FloatDPBounds(0.0_x, dp),
            FloatDPBounds(1.0_x, dp)});
        const Vector<FloatDPBounds> correlated_values = shared.field(same_input);
        const UpperBoxType correlated_image(correlated_values.size(), [&](SizeType j) {
            return UpperIntervalType(correlated_values[j].lower(), correlated_values[j].upper());
        });
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{1.0_x, 1.0_x}, {-1.0_x, -1.0_x}}),
            cast_exact_box(correlated_image))));

        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        using Model = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
        Scalar uncertain = Scalar::constant(
            domain, Interval<FloatDPUpperBound>(-1, 1, dp), sweeper);
        ValidatedIntervalTaylorBoxMultifunctionDP uncorrelated(Model({uncertain, uncertain}));
        const auto independent = uncorrelated.inclusion_integrator_field();
        ARIADNE_TESTING_EXPECT_EQUAL(independent.field.result_size(), 2u);
        ARIADNE_TESTING_PROVE(independent.inputs.dimension() >= 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(independent.field.argument_size(),
                                     2u + independent.inputs.dimension());
        const UpperBoxType state({{0.0_x, 0.0_x}, {0.0_x, 0.0_x}});
        const auto original = uncorrelated.evaluate(state);
        const Vector<FloatDPBounds> arguments(
            2u + independent.inputs.dimension(), [&](SizeType j) {
                if (j < 2u) {
                    return FloatDPBounds(0.0_x, dp);
                }
                const auto& input = independent.inputs[j - 2u];
                return FloatDPBounds(input.lower_bound(), input.upper_bound());
            });
        const auto expanded_values = independent.field(arguments);
        const UpperBoxType expanded(expanded_values.size(), [&](SizeType j) {
            return UpperIntervalType(expanded_values[j].lower(), expanded_values[j].upper());
        });
        ARIADNE_TESTING_PROVE(definitely(subset(original, cast_exact_box(expanded))));
        // (1,1) is allowed by the independent interval-coefficient box,
        // whereas the shared-input representation (u,-u) forbids it.
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{1.0_x, 1.0_x}, {1.0_x, 1.0_x}}), cast_exact_box(expanded))));
    }

    Void test_box_only_capability_does_not_select_solver() {
        // Both numerical solvers receive exactly the same box field.
        // The high-order solver must reject missing structured data, not
        // guess a model or silently dispatch to a different algorithm.
        class BoxOnlyAdapter : public ValidatedBoxMultifunctionInterface {
            ValidatedInputFieldBoxAdapterDP _wrapped;
          public:
            BoxOnlyAdapter(const EffectiveVectorMultivariateFunction& f,
                           const BoxDomainType& U, const BoxDomainType& domain)
                : _wrapped(f, U, domain) {}
            SizeType argument_size() const override { return _wrapped.argument_size(); }
            SizeType result_size() const override { return _wrapped.result_size(); }
            BoxDomainType domain() const override { return _wrapped.domain(); }
            UpperBoxType evaluate(const UpperBoxType& box) const override {
                return _wrapped.evaluate(box);
            }
        };

        BoxOnlyAdapter field_view(field(), inputs(), BoxDomainType({{-1.0_x, 1.0_x}}));
        ARIADNE_TESTING_PROVE(
            dynamic_cast<const InclusionIntegratorFieldProvider*>(&field_view) == nullptr);
        const ExactBoxType initial({{0.0_x, 0.125_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);
        FirstOrderMultifunctionIntegrator first_order;
        const MultifunctionIntegratorInterface& low = first_order;
        const auto boxes = low.flow_step(field_view, initial, h, bound);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(boxes.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(boxes.evolve)));

        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator high_order(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const MultifunctionIntegratorInterface& high = high_order;
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            high.flow_step(field_view, initial, h, bound));
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
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_every_endpoint_strategy_encloses_full_input_range);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_invalid_time_step);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_parametric_box_extension_and_low_order_bridge);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_high_order_common_integrator_bridge);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_high_order_interval_taylor_via_parametric_conversion);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_interval_taylor_conversion_of_state_dependent_input);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_interval_taylor_conversion_preserves_switching_endpoints);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_multidimensional_interval_taylor_conversion);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_shared_input_correlation_is_explicit);
    ARIADNE_TESTING_ADD_TEST(TestSelectedInclusionIntegrator, test_box_only_capability_does_not_select_solver);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
