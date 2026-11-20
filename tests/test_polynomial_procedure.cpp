/***************************************************************************
 *            test_polynomial_procedure.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file tests/test_polynomial_procedure.cpp
 *  @brief Test the private polynomial recurrence via its public integrator.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#include "function/function.hpp"
#include "integration/integrator.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

// Capture the integrator's existing diagnostics without exporting its
// implementation helpers for the sole purpose of testing.
template<class F>
auto capture_diagnostics(F&& operation) {
    std::ostringstream log;
    std::streambuf* const original = std::cerr.rdbuf(log.rdbuf());
    try {
        auto result = operation();
        std::cerr.rdbuf(original);
        return std::make_pair(std::move(result), log.str());
    } catch (...) {
        std::cerr.rdbuf(original);
        throw;
    }
}

auto integrate_with_diagnostics(
    const ValidatedVectorMultivariateFunction& field,
    const ExactBoxType& domain) {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error = 1e-6);
    integrator.set_preconditioning(TaylorSeriesPreconditioning::QR);
    integrator.set_minimum_spacial_order(2u);
    integrator.set_maximum_spacial_order(2u);
    integrator.set_minimum_temporal_order(2u);
    integrator.set_maximum_temporal_order(2u);
    integrator.set_diagnostics(true);
    return capture_diagnostics([&] {
        const auto initial = integrator.function_factory().create_identity(domain);
        return integrator.step(field, integrator.precondition(initial),
                               suggest(StepSizeType(0.00390625_x)));
    });
}

Bool diagnostic_reports(const std::string& log, const std::string& value) {
    const auto start = log.find("[CentrePolynomialDefectDiagnostic]");
    if (start == std::string::npos) {
        return false;
    }
    const auto end = log.find('\n', start);
    return log.substr(start, end - start).find(value) != std::string::npos;
}

class TestPolynomialProcedure : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_polynomial_operators() {
        // Polynomial including negation, subtraction, halving/squaring and
        // division by a constant. Exercise the same QR recurrence that
        // invokes polynomial_procedure_degree and Differential evaluation.
        const auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const ValidatedVectorMultivariateFunction field =
            EffectiveVectorMultivariateFunction({-x + hlf(sqr(x)) - x / 4});
        const ExactBoxType domain({{1.875_x, 2.125_x}});
        const auto [step, log] = integrate_with_diagnostics(field, domain);
        ARIADNE_TESTING_PROVE(
            diagnostic_reports(log, "exact_polynomial_available=1"));
        ARIADNE_TESTING_PROVE(
            diagnostic_reports(log, "exact_polynomial_degree="));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.flowpipe_mapping().range())));
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolved_mapping().range())));

        // x'=x(x-5/2)/2, x(0)=2:
        // x(t)=(5/2)/(1+(1/4)*exp((5/4)*t)).
        const FloatDPBounds h(step.time_step(), dp);
        const FloatDPBounds expected =
            FloatDPBounds(2.5_x, dp) /
            (FloatDPBounds(1, dp) + FloatDPBounds(0.25_x, dp) *
             exp(FloatDPBounds(1.25_x, dp) * h));
        const Vector<FloatDPBounds> point({FloatDPBounds(2, dp)});
        ARIADNE_TESTING_PROVE(
            refines(expected, step.evolved_mapping()(point)[0u]));
    }

    Void test_unsupported_polynomial_diagnostics() {
        // These operators are not accepted by the internal polynomial
        // analyser. The solver must nevertheless return rigorous flows.
        const auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const ValidatedVectorMultivariateFunction reciprocal =
            EffectiveVectorMultivariateFunction({1 / x});
        const ValidatedVectorMultivariateFunction cubic =
            EffectiveVectorMultivariateFunction({pow(x, 3)});
        const ValidatedVectorMultivariateFunction sine =
            EffectiveVectorMultivariateFunction({sin(x)});
        const ExactBoxType domain({{0.875_x, 1.125_x}});
        for (const auto& field: {reciprocal, cubic, sine}) {
            const auto [step, log] = integrate_with_diagnostics(field, domain);
            ARIADNE_TESTING_PROVE(
                diagnostic_reports(log, "exact_polynomial_available=0"));
            ARIADNE_TESTING_PROVE(definitely(is_bounded(step.flowpipe_mapping().range())));
            ARIADNE_TESTING_PROVE(definitely(is_bounded(step.evolved_mapping().range())));
        }
    }

    Void test_non_polynomial_recurrence_diagnostics() {
        const auto x = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
        const ValidatedVectorMultivariateFunction field =
            EffectiveVectorMultivariateFunction({sin(x)});
        const ExactBoxType domain({{-0.125_x, 0.125_x}});
        const auto [step, log] = integrate_with_diagnostics(field, domain);
        ARIADNE_TESTING_PROVE(
            diagnostic_reports(log, "exact_polynomial_available=0"));
        ARIADNE_TESTING_PROVE(
            diagnostic_reports(log, "recurrence_defect_range="));
        ARIADNE_TESTING_EXPECT_EQUAL(step.evolved_mapping().argument_size(), 1u);
        ARIADNE_TESTING_PROVE(definitely(is_bounded(step.flowpipe_mapping().range())));
        // x=0 is an equilibrium of the nonpolynomial ODE.
        const Vector<FloatDPBounds> point({FloatDPBounds(0, dp)});
        ARIADNE_TESTING_PROVE(
            refines(FloatDPBounds(0, dp), step.evolved_mapping()(point)[0u]));
    }
};

} // namespace

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_polynomial_procedure");
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_polynomial_operators);
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_unsupported_polynomial_diagnostics);
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_non_polynomial_recurrence_diagnostics);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
