/***************************************************************************
 *            test_multifunction_integrator.cpp
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

/** @file tests/test_multifunction_integrator.cpp
 *  @brief Tests for validated box integration of differential inclusions.
 */

#include <cstdlib>
#include <stdexcept>
#include <string>

#include "function/taylor_multifunction.hpp"
#include "integration/multifunction_integrator.hpp"

#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

class TestMultifunctionIntegrator : public Ariadne::Testing::ContextFixture {
    using Scalar = ValidatedIntervalTaylorFunctionModel<FloatDP>;
    using VectorModel = ValidatedVectorIntervalTaylorFunctionModel<FloatDP>;
    using Ivl = Interval<FloatDPUpperBound>;

  public:
    using ContextFixture::ContextFixture;

    Void test_constant_inclusion() {
        BoxDomainType domain({{-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar uncertain = Scalar::constant(domain, Ivl(-1, 1, dp), sweeper);
        VectorModel model({uncertain});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        FirstOrderMultifunctionIntegrator integrator;
        const MultifunctionIntegratorInterface& polymorphic = integrator;

        ExactBoxType initial({{0, 0}});
        UpperBoxType bound({{-1, 1}});
        StepSizeType h(0.25_x);
        const auto step = polymorphic.flow_step(field, initial, h, bound);
        ARIADNE_TESTING_EXPECT_EQUAL(step.time_step, h);
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{-0.25_x, 0.25_x}}), step.reach)));
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{-0.25_x, 0.25_x}}), step.evolve)));
        ARIADNE_TESTING_PROVE(definitely(inside(step.reach, bound)));
    }

    Void test_time_varying_selection() {
        // x' = u(t), y' = x, u(t) in [-1,1]. A switching input with
        // u=+1 then -1 produces (x(1/2),y(1/2))=(0,1/16).
        BoxDomainType domain({{-1, 1}, {-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar x = Scalar::coordinate(domain, 0u, sweeper);
        Scalar uncertain = Scalar::constant(domain, Ivl(-1, 1, dp), sweeper);
        VectorModel model({uncertain, x});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        FirstOrderMultifunctionIntegrator integrator;
        ExactBoxType initial({{0, 0}, {0, 0}});
        UpperBoxType bound({{-1, 1}, {-1, 1}});

        const auto step = integrator.flow_step(field, initial, StepSizeType(0.5_x), bound);
        const ExactBoxType switching_endpoint({{0, 0}, {0.0625_x, 0.0625_x}});
        ARIADNE_TESTING_PROVE(definitely(subset(switching_endpoint, step.evolve)));
        ARIADNE_TESTING_PROVE(definitely(subset(switching_endpoint, step.reach)));
    }

    Void test_reject_uncertified_bound() {
        BoxDomainType domain({{-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar uncertain = Scalar::constant(domain, Ivl(-1, 1, dp), sweeper);
        VectorModel model({uncertain});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        FirstOrderMultifunctionIntegrator integrator;
        ExactBoxType initial({{0, 0}});

        // The reach candidate cannot fit inside this state bound.
        UpperBoxType insufficient({{-0.125_x, 0.125_x}});
        ARIADNE_TESTING_EXPECT_THROWS(FlowBoundsException, integrator.flow_step(field, initial, StepSizeType(0.5_x), insufficient));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, integrator.flow_step(field, initial, StepSizeType(0u), UpperBoxType({{-1, 1}})));

        // A box extension must not be evaluated beyond the model domain.
        UpperBoxType outside_domain({{-2, 2}});
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, field.evaluate(outside_domain));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, integrator.flow_step(field, initial, StepSizeType(0.5_x), outside_domain));
    }

    // An interior bound must certify the very same constant inclusion.
    // This isolates the continuation test from domain-boundary handling.
    Void test_interior_valid_step() {
        BoxDomainType domain({{-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar uncertain = Scalar::constant(domain, Ivl(-1, 1, dp), sweeper);
        VectorModel model({uncertain});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        FirstOrderMultifunctionIntegrator integrator;
        ExactBoxType initial({{0, 0}});
        UpperBoxType interior({{-0.5_x, 0.5_x}});
        const auto step = integrator.flow_step(field, initial, StepSizeType(0.25_x), interior);
        ARIADNE_TESTING_PROVE(definitely(inside(step.reach, interior)));
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{-0.25_x, 0.25_x}}), step.evolve)));
    }

    // The equality case is legitimate: the bounding box may coincide with
    // the (closed) domain of the Taylor model.
    Void test_boundary_domain_check() {
        BoxDomainType domain({{-1, 1}});
        UpperBoxType same_bound({{-1, 1}});
        ARIADNE_TESTING_PROVE(definitely(subset(same_bound, domain)));
    }

    Void test_boundary_box_evaluation() {
        BoxDomainType domain({{-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar uncertain = Scalar::constant(domain, Ivl(-1, 1, dp), sweeper);
        VectorModel model({uncertain});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        UpperBoxType same_bound({{-1, 1}});
        const auto values = field.evaluate(same_bound);
        ARIADNE_TESTING_EXPECT_EQUAL(values.dimension(), 1u);
        ARIADNE_TESTING_PROVE(definitely(subset(ExactBoxType({{-1, 1}}), values)));
    }

    Void test_reject_dimension_mismatch() {
        BoxDomainType domain({{-1, 1}});
        ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        Scalar x = Scalar::coordinate(domain, 0u, sweeper);
        VectorModel model({x, x});
        ValidatedIntervalTaylorBoxMultifunctionDP field(model);
        FirstOrderMultifunctionIntegrator integrator;
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, integrator.flow_step(field, ExactBoxType({{0, 0}}), StepSizeType(0.25_x), UpperBoxType({{-1, 1}})));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_multifunction_integrator");
    // Diagnostic selection is opt-in; production runs exercise all four cases.
    const char* requested = std::getenv("ARIADNE_MULTIFUNCTION_CASE");
    const std::string selected = requested ? requested : "";
    const auto selected_case = [&](const char* key) {
        return selected.empty() || ("," + selected + ",").find("," + std::string(key) + ",") != std::string::npos;
    };
    if (selected_case("constant_inclusion")) {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_constant_inclusion);
    }
    if (selected_case("time_varying_selection")) {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_time_varying_selection);
    }
    if (selected_case("reject_uncertified_bound")) {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_reject_uncertified_bound);
    }
    if (selected_case("reject_dimension_mismatch")) {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_reject_dimension_mismatch);
    }
    if (selected.empty() || selected == "interior_valid_step") {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_interior_valid_step);
    }
    if (selected.empty() || selected == "boundary_domain_check") {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_boundary_domain_check);
    }
    if (selected.empty() || selected == "boundary_box_evaluation") {
        ARIADNE_TESTING_ADD_TEST(TestMultifunctionIntegrator, test_boundary_box_evaluation);
    }
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
