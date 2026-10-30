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

namespace Ariadne::IntegrationDetail {
// Internal polynomial diagnostics, shared by the recurrence and its unit tests.
// DIV requires a constant denominator; POW and non-polynomial operations are
// deliberately excluded from the degree analyser.
Bool polynomial_procedure_degree(Vector<ValidatedProcedure> const& p, DegreeType& degree);
Bool evaluate_polynomial_procedure(Vector<ValidatedProcedure> const& p,
    Vector<Differential<FloatDPBounds>> const& x,
    Vector<Differential<FloatDPBounds>>& result);
}

#endif
