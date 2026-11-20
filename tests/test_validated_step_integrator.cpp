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
        const ValidatedStepIntegratorInterface& solvers[2] = {flow_solver, inclusion_solver};
        const auto ordinary_step = solvers[0].step(function_request);
        const auto inclusion_step = solvers[1].step(multifunction_request);

        ARIADNE_TESTING_EXPECT_EQUAL(ordinary_step.time_step(), h);
        ARIADNE_TESTING_EXPECT_EQUAL(inclusion_step.time_step(), h);
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{0.03125_x, 0.15625_x}}),
            cast_exact_box(ordinary_step.evolve()))));
        ARIADNE_TESTING_PROVE(definitely(subset(
            ExactBoxType({{-0.03125_x, 0.15625_x}}),
            cast_exact_box(inclusion_step.evolve()))));
        ARIADNE_TESTING_PROVE(
            ordinary_step.representation_as<FlowStepModelType>() != nullptr);
        ARIADNE_TESTING_PROVE(not inclusion_step.has_representation());

        // An unsuitable numerical input must not select another solver.
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            solvers[0].step(multifunction_request));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            solvers[1].step(function_request));
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
        const auto* original = step.representation_as<SelectedInclusionStep>();
        ARIADNE_TESTING_PROVE(original != nullptr);
        ARIADNE_TESTING_PROVE(original->reach_models.size() >= 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(original->time_step, h);
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
    ARIADNE_TESTING_ADD_TEST(TestValidatedStepIntegrator, test_reject_invalid_request);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
