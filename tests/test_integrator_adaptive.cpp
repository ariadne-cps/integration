/***************************************************************************
 *            test_integrator.cpp
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

#include <string>
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
Void test_fixed_step_partial_flow() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({x});
    TaylorSeriesIntegrator integrator(Sweeper<FloatDP>(),lipschitz_tolerance=0.5_x,order=2);
    ExactBoxType domain={{1.0_x,1.0_x}};
    StepSizeType h(1);
    try {
        integrator.IntegratorBase::flow_step(field,domain,h);
        ARIADNE_TEST_ERROR("Expected a partial flow after reducing the requested fixed step")
    } catch(IncompleteFlowException const& error) {
        auto const& flow=error.computed_model();
        auto time=flow.domain()[1u].upper_bound();
        ARIADNE_TEST_COMPARE(time,<,FloatDP(h,dp))
        ARIADNE_TEST_COMPARE(time,>,FloatDP(0,dp))
        Vector<FloatDPBounds> point({FloatDPBounds(1,dp),FloatDPBounds(time)});
        ARIADNE_TEST_ASSERT(refines(exp(FloatDPBounds(time)),flow(point)[0u]))
    }
}

Void test_fixed_step_reduction_limit() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    // A bounded but rapidly varying field needs more than 32 halvings for
    // this low-order enclosure to fit the fixed a-priori bounding box.
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({sin(0x1p100_x*x)});
    TaylorSeriesIntegrator integrator(Sweeper<FloatDP>(),lipschitz_tolerance=0.5_x,order=2);
    ExactBoxType domain={{1.0_x,1.0_x}};
    try {
        integrator.IntegratorBase::flow_step(field,domain,StepSizeType(1));
        ARIADNE_TEST_ERROR("Expected the fixed-step reduction limit")
    } catch(IncompleteFlowException const& error) {
        ARIADNE_TEST_ASSERT(std::string(error.what()).find("32 step halvings")!=std::string::npos)
        ARIADNE_TEST_EQUALS(error.computed_model().domain()[1u].upper_bound(),FloatDP(0x1p-32_x,dp))
    }
}

Void test_suggested_base_step_reduction() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({-x});
    ExactBoxType domain={{-1.0_x,1.0_x}};
    TaylorSeriesIntegrator integrator(step_maximum_error=1e-6,order=2);
    auto flow=integrator.IntegratorBase::flow_step(field,domain,suggest(StepSizeType(1)));
    auto time=flow.domain()[1u].upper_bound();
    ARIADNE_TEST_COMPARE(time,<,FloatDP(1,dp))
    ARIADNE_TEST_COMPARE(time,>,FloatDP(0,dp))
    ARIADNE_TEST_COMPARE(flow.error(),<=,1e-6_pr)
    Vector<FloatDPBounds> point({FloatDPBounds(1,dp),FloatDPBounds(time)});
    ARIADNE_TEST_ASSERT(refines(exp(-FloatDPBounds(time)),flow(point)[0u]))
}
}

int main() {
    ARIADNE_TEST_CALL(test_fixed_step_partial_flow())
    ARIADNE_TEST_CALL(test_fixed_step_reduction_limit())
    ARIADNE_TEST_CALL(test_suggested_base_step_reduction())
    return ARIADNE_TEST_FAILURES;
}

