/***************************************************************************
 *            test_polynomial_procedure.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "polynomial_procedure.hpp"
#include "function/formula.hpp"
#include "numeric/numeric.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;
using namespace Ariadne::IntegrationDetail;

namespace {
Void test_polynomial_operators() {
    auto x=ValidatedFormula::coordinate(0u);
    auto positive=ValidatedFormula::unary(Pos(),x);
    auto negative=ValidatedFormula::unary(Neg(),positive);
    auto half_square=ValidatedFormula::unary(Hlf(),ValidatedFormula::unary(Sqr(),x));
    auto quarter=ValidatedFormula::binary(Div(),x,ValidatedFormula::constant(ValidatedNumber(4)));
    auto expression=ValidatedFormula::binary(Sub(),ValidatedFormula::binary(Add(),negative,half_square),quarter);
    Vector<ValidatedProcedure> procedure(1u,Vector<ValidatedFormula>({expression}));
    DegreeType degree=0u;
    ARIADNE_TEST_ASSERT(polynomial_procedure_degree(procedure,degree))
    ARIADNE_TEST_EQUALS(degree,2u)
    using D=Differential<FloatDPBounds>;
    Vector<D> arguments({D::variable(1u,3u,FloatDPBounds(2,dp),0u)});
    Vector<D> evaluated(1u,arguments.zero_element());
    ARIADNE_TEST_ASSERT(evaluate_polynomial_procedure(procedure,arguments,evaluated))
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(-0.5_x,dp),evaluated[0u].value()))
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(0.75_x,dp),evaluated[0u][MultiIndex::unit(1u,0u)]))
    MultiIndex second(1u);second[0u]=2u;
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(0.5_x,dp),evaluated[0u][second]))
}

Void test_unsupported_polynomial_diagnostics() {
    auto x=ValidatedFormula::coordinate(0u);
    auto variable_denominator=ValidatedFormula::binary(Div(),ValidatedFormula::constant(ValidatedNumber(1)),x);
    using D=Differential<FloatDPBounds>;
    Vector<D> arguments({D::variable(1u,3u,FloatDPBounds(2,dp),0u)});
    Vector<D> evaluated(1u,arguments.zero_element());
    for(auto expression:{variable_denominator,pow(x,3),sin(x)}) {
        Vector<ValidatedProcedure> procedure(1u,Vector<ValidatedFormula>({expression}));
        DegreeType degree=0u;
        ARIADNE_TEST_ASSERT(not polynomial_procedure_degree(procedure,degree))
    }
    Vector<ValidatedProcedure> non_polynomial(1u,Vector<ValidatedFormula>({sin(x)}));
    ARIADNE_TEST_ASSERT(not evaluate_polynomial_procedure(non_polynomial,arguments,evaluated))
}
}

int main() {
    ARIADNE_TEST_CALL(test_polynomial_operators())
    ARIADNE_TEST_CALL(test_unsupported_polynomial_diagnostics())
    return ARIADNE_TEST_FAILURES;
}

