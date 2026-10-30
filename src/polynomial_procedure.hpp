/***************************************************************************
 *            polynomial_procedure.hpp
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

#ifndef ARIADNE_INTEGRATION_POLYNOMIAL_PROCEDURE_HPP
#define ARIADNE_INTEGRATION_POLYNOMIAL_PROCEDURE_HPP

#include "numeric/numeric.hpp"
#include "algebra/differential.hpp"
#include "function/procedure.hpp"
#include "function/taylor_function.hpp"

namespace Ariadne::IntegrationDetail {
// Internal polynomial diagnostics, shared by the recurrence and its unit tests.
// DIV requires a constant denominator; POW and non-polynomial operations are
// deliberately excluded from the degree analyser.
Bool polynomial_procedure_degree(Vector<ValidatedProcedure> const& p, DegreeType& degree);
Bool evaluate_polynomial_procedure(Vector<ValidatedProcedure> const& p,
    Vector<Differential<FloatDPBounds>> const& x,
    Vector<Differential<FloatDPBounds>>& result);

// A centre-series candidate and residual diagnostics. The candidate receives
// a validated remainder only in the caller after separate certification.
struct CentrePolynomialRecurrenceResult {
    ValidatedVectorMultivariateTaylorFunctionModelDP polynomial;
    ValidatedVectorMultivariateTaylorFunctionModelDP defect;
    Vector<FloatDPBounds> direct_defect_range;
    Bool exact_polynomial_available;
    DegreeType exact_polynomial_degree;
    Vector<FloatDPBounds> exact_polynomial_defect_range;
};
CentrePolynomialRecurrenceResult graded_series_centre_polynomial_step(
    Vector<ValidatedProcedure> const& p, ExactBoxType const& domain,
    ExactIntervalType const& time, Sweeper<FloatDP> const& sweeper,
    DegreeType spatial_order, DegreeType temporal_order, Bool exact_diagnostics);

}

#endif
