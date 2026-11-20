/***************************************************************************
 *            test_declared_multifunction.cpp
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

/** @file test_declared_multifunction.cpp
 *  @brief Test explicit constant-parameter and time-input semantics.
 */

#include <stdexcept>

#include "function/formula.tpl.hpp"
#include "integration/declared_multifunction.hpp"
#include "integration/selected_inclusion_integrator.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestDeclaredMultifunction : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    EffectiveVectorMultivariateFunction mixed_field() {
        // One state coordinate, one fixed parameter and one time input.
        auto p = EffectiveFormula::coordinate(1u);
        auto u = EffectiveFormula::coordinate(2u);
        return EffectiveVectorMultivariateFormulaFunction(
            3u, Vector<EffectiveFormula>{p + u});
    }

    BoxDomainType states() { return BoxDomainType({{-1.0_x, 1.0_x}}); }
    BoxDomainType parameters() { return BoxDomainType({{-0.5_x, 0.5_x}}); }
    BoxDomainType inputs() { return BoxDomainType({{-0.25_x, 0.25_x}}); }

    Void test_distinct_roles_and_first_order_inclusion() {
        ValidatedDeclaredBoxMultifunctionDP declared(mixed_field(), states(), parameters(), inputs());
        ARIADNE_TESTING_EXPECT_EQUAL(declared.argument_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(declared.result_size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(declared.constant_parameters().dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(declared.time_inputs().dimension(), 1u);
        ARIADNE_TESTING_PROVE(not declared.parametric_form().has_value());

        const auto image = declared.evaluate(UpperBoxType({{0.0_x, 0.0_x}}));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.75_x, 0.75_x}}), cast_exact_box(image))));
        const FirstOrderMultifunctionIntegrator integrator;
        const auto step = integrator.flow_step(
            declared, ExactBoxType({{0.0_x, 0.0_x}}),
            StepSizeType(0.03125_x), UpperBoxType({{-1.0_x, 1.0_x}}));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.0234375_x, 0.0234375_x}}), cast_exact_box(step.evolve))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0_x, 0.0_x}}), cast_exact_box(step.reach))));
    }

    Void test_high_order_requires_explicit_relaxation() {
        ValidatedDeclaredBoxMultifunctionDP declared(mixed_field(), states(), parameters(), inputs());
        TaylorSeriesIntegrator base(ThresholdSweeper<FloatDP>(dp, 1e-10),
                                    lipschitz_tolerance = 0.5_x, order = 8);
        SelectedParametricMultifunctionIntegrator selected(
            base, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation()));
        const ExactBoxType initial({{0.0_x, 0.0_x}});
        const UpperBoxType bound({{-1.0_x, 1.0_x}});
        const StepSizeType h(0.03125_x);

        // A frozen trajectory parameter is never silently reclassified.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            selected.flow_step(declared, initial, h, bound));

        // An explicit over-approximation may be chosen for this backend.
        const auto relaxation = declared.time_varying_relaxation();
        ARIADNE_TESTING_EXPECT_EQUAL(relaxation.inputs.dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(relaxation.field.argument_size(), 3u);
        ValidatedParametricBoxMultifunctionDP relaxed_field(
            relaxation.field, relaxation.inputs, declared.domain());
        const auto step = selected.flow_step(relaxed_field, initial, h, bound);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.reach)));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolve)));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.0_x, 0.0_x}}), cast_exact_box(step.reach))));
    }

    Void test_only_time_inputs_are_directly_integrable() {
        auto u = EffectiveFormula::coordinate(1u);
        EffectiveVectorMultivariateFunction f = EffectiveVectorMultivariateFormulaFunction(
            2u, Vector<EffectiveFormula>{u});
        ValidatedDeclaredBoxMultifunctionDP declared(f, states(), BoxDomainType(0u), inputs());
        const auto form = declared.parametric_form();
        ARIADNE_TESTING_PROVE(form.has_value());
        ARIADNE_TESTING_EXPECT_EQUAL(form->field.argument_size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(form->inputs.dimension(), 1u);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.25_x, 0.25_x}}),
            cast_exact_box(declared.evaluate(UpperBoxType({{0.0_x, 0.0_x}}))))));
    }

    Void test_shared_stationary_parameter_is_not_split() {
        BoxDomainType state_domain({{-1.0_x, 1.0_x}, {-1.0_x, 1.0_x}});
        auto p = EffectiveFormula::coordinate(2u);
        EffectiveVectorMultivariateFunction f = EffectiveVectorMultivariateFormulaFunction(
            3u, Vector<EffectiveFormula>{p, -p});
        ValidatedDeclaredBoxMultifunctionDP declared(
            f, state_domain, BoxDomainType({{-1.0_x, 1.0_x}}), BoxDomainType(0u));
        ARIADNE_TESTING_PROVE(not declared.parametric_form().has_value());
        const auto relaxed = declared.time_varying_relaxation();
        ARIADNE_TESTING_EXPECT_EQUAL(relaxed.inputs.dimension(), 1u);

        const Vector<FloatDPBounds> arguments({
            FloatDPBounds(0.0_x, dp),
            FloatDPBounds(0.0_x, dp),
            FloatDPBounds(1.0_x, dp)});
        const auto values = relaxed.field(arguments);
        ARIADNE_TESTING_EXPECT_EQUAL(values.size(), 2u);
        const UpperBoxType output(values.size(), [&](SizeType i) {
            return UpperIntervalType(values[i].lower(), values[i].upper());
        });
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{1.0_x, 1.0_x}, {-1.0_x, -1.0_x}}),
            cast_exact_box(output))));
    }

    Void test_invalid_dimensions_and_parameter_free_ode() {
        auto x = EffectiveFormula::coordinate(0u);
        EffectiveVectorMultivariateFunction ode = EffectiveVectorMultivariateFormulaFunction(
            1u, Vector<EffectiveFormula>{x});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedDeclaredBoxMultifunctionDP(ode, states(), parameters(), inputs()));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            ValidatedDeclaredBoxMultifunctionDP(ode, BoxDomainType(0u),
                                                BoxDomainType(0u), BoxDomainType(0u)));
        ValidatedDeclaredBoxMultifunctionDP declared(
            ode, states(), BoxDomainType(0u), BoxDomainType(0u));
        ARIADNE_TESTING_PROVE(not declared.parametric_form().has_value());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            declared.time_varying_relaxation());
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            declared.evaluate(UpperBoxType({{-2.0_x, 2.0_x}})));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_declared_multifunction");
    ARIADNE_TESTING_ADD_TEST(TestDeclaredMultifunction, test_distinct_roles_and_first_order_inclusion);
    ARIADNE_TESTING_ADD_TEST(TestDeclaredMultifunction, test_high_order_requires_explicit_relaxation);
    ARIADNE_TESTING_ADD_TEST(TestDeclaredMultifunction, test_only_time_inputs_are_directly_integrable);
    ARIADNE_TESTING_ADD_TEST(TestDeclaredMultifunction, test_shared_stationary_parameter_is_not_split);
    ARIADNE_TESTING_ADD_TEST(TestDeclaredMultifunction, test_invalid_dimensions_and_parameter_free_ode);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
