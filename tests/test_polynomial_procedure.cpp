/***************************************************************************
 *            test_polynomial_procedure.cpp
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

/** @file tests/test_polynomial_procedure.cpp
 *  @brief Test polynomial procedure behavior.
 */

#include "polynomial_procedure.hpp"
#include "function/formula.hpp"
#include "numeric/numeric.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace Ariadne::IntegrationDetail;

namespace {
class TestPolynomialProcedure : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

Void test_polynomial_operators() {
    auto x=ValidatedFormula::coordinate(0u);
    auto positive=ValidatedFormula::unary(Pos(),x);
    auto negative=ValidatedFormula::unary(Neg(),positive);
    auto half_square=ValidatedFormula::unary(Hlf(),ValidatedFormula::unary(Sqr(),x));
    auto quarter=ValidatedFormula::binary(Div(),x,ValidatedFormula::constant(ValidatedNumber(4)));
    auto expression=ValidatedFormula::binary(Sub(),ValidatedFormula::binary(Add(),negative,half_square),quarter);
    Vector<ValidatedProcedure> procedure(1u,Vector<ValidatedFormula>({expression}));
    DegreeType polynomial_degree=0u;
    ARIADNE_TESTING_PROVE(polynomial_procedure_degree(procedure,polynomial_degree));
    ARIADNE_TESTING_EXPECT_EQUAL(polynomial_degree,2u);
    using D=Differential<FloatDPBounds>;
    Vector<D> arguments({D::variable(1u,3u,FloatDPBounds(2,dp),0u)});
    Vector<D> evaluated(1u,arguments.zero_element());
    ARIADNE_TESTING_PROVE(evaluate_polynomial_procedure(procedure,arguments,evaluated));
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(-0.5_x,dp),evaluated[0u].value()));
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.75_x,dp),evaluated[0u][MultiIndex::unit(1u,0u)]));
    MultiIndex second(1u);second[0u]=2u;
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.5_x,dp),evaluated[0u][second]));
}

Void test_unsupported_polynomial_diagnostics() {
    auto x=ValidatedFormula::coordinate(0u);
    auto variable_denominator=ValidatedFormula::binary(Div(),ValidatedFormula::constant(ValidatedNumber(1)),x);
    using D=Differential<FloatDPBounds>;
    Vector<D> arguments({D::variable(1u,3u,FloatDPBounds(2,dp),0u)});
    Vector<D> evaluated(1u,arguments.zero_element());
    for(auto expression:{variable_denominator,pow(x,3),sin(x)}) {
        Vector<ValidatedProcedure> procedure(1u,Vector<ValidatedFormula>({expression}));
        DegreeType polynomial_degree=0u;
        ARIADNE_TESTING_PROVE(not polynomial_procedure_degree(procedure,polynomial_degree));
    }
    Vector<ValidatedProcedure> non_polynomial(1u,Vector<ValidatedFormula>({sin(x)}));
    ARIADNE_TESTING_PROVE(not evaluate_polynomial_procedure(non_polynomial,arguments,evaluated));
}

Void test_non_polynomial_recurrence_diagnostics() {
    auto x=ValidatedFormula::coordinate(0u);
    Vector<ValidatedProcedure> field(1u,Vector<ValidatedFormula>({sin(x)}));
    ExactBoxType domain={{-0.125_x,0.125_x}};
    ExactIntervalType time(0.0_x,0.03125_x);
    auto recurrence=graded_series_centre_polynomial_step(field,domain,time,Sweeper<FloatDP>(),2u,2u,true);
    ARIADNE_TESTING_PROVE(not recurrence.exact_polynomial_available);
    ARIADNE_TESTING_EXPECT_EQUAL(recurrence.polynomial.result_size(),1u);
    ARIADNE_TESTING_EXPECT_EQUAL(recurrence.polynomial.argument_size(),2u);
    Vector<FloatDPBounds> point({FloatDPBounds(0,dp),FloatDPBounds(0.03125_x,dp)});
    // sin(0)=0: the equilibrium has zero candidate value and zero residual.
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0,dp),recurrence.polynomial(point)[0u]));
    ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0,dp),recurrence.defect(point)[0u]));
    ARIADNE_TESTING_PROVE(definitely(is_bounded(recurrence.defect.range())));
}
};

}


int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_polynomial_procedure");
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_polynomial_operators);
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_unsupported_polynomial_diagnostics);
    ARIADNE_TESTING_ADD_TEST(TestPolynomialProcedure, test_non_polynomial_recurrence_diagnostics);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
