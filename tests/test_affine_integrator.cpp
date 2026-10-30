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

#include <memory>
#include <sstream>
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
Void test_linear_flow() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({x});
    ExactBoxType domain={{0.125_x,0.375_x}};
    UpperBoxType bound={{0.0_x,0.5_x}};
    StepSizeType h(0.125_x);
    AffineIntegrator integrator(SpacialOrder(1u),TemporalOrder(6u));
    ARIADNE_TEST_EQUALS(integrator.spacial_order(),1u)
    ARIADNE_TEST_EQUALS(integrator.temporal_order(),6u)
    std::unique_ptr<IntegratorInterface> clone(integrator.clone());
    std::ostringstream text;
    text << *clone;
    ARIADNE_TEST_ASSERT(text.str().find("AffineIntegrator")!=std::string::npos)
    auto flow=clone->flow_step(field,domain,h,bound);
    Vector<FloatDPBounds> point({FloatDPBounds(0.25_x,dp),FloatDPBounds(h,dp)});
    auto expected=FloatDPBounds(0.25_x,dp)*exp(FloatDPBounds(h,dp));
    ARIADNE_TEST_ASSERT(refines(expected,flow(point)[0u]))
}

Void test_nonlinear_flow_enclosure() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({sqr(x)});
    ExactBoxType domain={{0.125_x,0.375_x}};
    UpperBoxType bound={{0.0_x,0.5_x}};
    StepSizeType h(0.125_x);
    AffineIntegrator integrator(SpacialOrder(1u),TemporalOrder(6u));
    auto flow=integrator.flow_step(field,domain,h,bound);
    // x'=x^2 has solution x0/(1-x0*t). Check initial, interior and
    // endpoint times at both spatial endpoints and at the centre.
    for(ExactDouble initial:{0.125_x,0.25_x,0.375_x}) {
      for(ExactDouble time:{0.0_x,0.0625_x,0.125_x}) {
        FloatDPBounds x0(initial,dp),t(time,dp);
        Vector<FloatDPBounds> point({x0,t});
        auto expected=x0/(1-x0*t);
        ARIADNE_TEST_ASSERT(refines(expected,flow(point)[0u]))
      }
    }
    ARIADNE_TEST_ASSERT(flow.error().get(dp).raw()>FloatDP(0,dp))
}

Void test_time_and_parameter_flow() {
    auto t=EffectiveScalarMultivariateFunction::coordinate(3u,1u);
    auto a=EffectiveScalarMultivariateFunction::coordinate(3u,2u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({t+a});
    ExactBoxType domain={{0.125_x,0.375_x}};
    ExactBoxType parameters={{-0.125_x,0.125_x}};
    Interval<StepSizeType> time(StepSizeType(1),StepSizeType(1.125_x));
    UpperBoxType bound={{0.0_x,1.0_x}};
    AffineIntegrator integrator(SpacialOrder(1u),TemporalOrder(6u));
    auto flow=integrator.flow_step(field,domain,time,parameters,bound);
    for(ExactDouble initial:{0.125_x,0.375_x}) {
      for(ExactDouble parameter:{-0.125_x,0.125_x}) {
        Vector<FloatDPBounds> point({FloatDPBounds(initial,dp),FloatDPBounds(1.125_x,dp),FloatDPBounds(parameter,dp)});
        auto expected=FloatDPBounds(initial,dp)+FloatDPBounds(parameter,dp)*FloatDPBounds(0.125_x,dp)
                      +(sqr(FloatDPBounds(1.125_x,dp))-1)/2;
        ARIADNE_TEST_ASSERT(refines(expected,flow(point)[0u]))
      }
    }
}

Void test_flow_derivatives() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({sqr(x)});
    AffineIntegrator integrator(SpacialOrder(1u),TemporalOrder(3u));
    auto derivatives=integrator.flow_derivative(field,Vector<FloatDPBounds>({FloatDPBounds(0.25_x,dp)}));
    ARIADNE_TEST_EQUALS(derivatives.size(),1u)
    ARIADNE_TEST_EQUALS(derivatives[0u].argument_size(),2u)
    FloatDPBounds temporal_coefficient=derivatives[0u][MultiIndex::unit(2u,1u)];
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(0.0625_x,dp),temporal_coefficient))
    AffineIntegrator stationary(SpacialOrder(0u),TemporalOrder(0u));
    auto zero_order=stationary.flow_derivative(field,Vector<FloatDPBounds>({FloatDPBounds(0.25_x,dp)}));
    ARIADNE_TEST_ASSERT(refines(FloatDPBounds(0.25_x,dp),zero_order[0u].value()))
}
}

int main() {
    ARIADNE_TEST_CALL(test_linear_flow())
    ARIADNE_TEST_CALL(test_nonlinear_flow_enclosure())
    ARIADNE_TEST_CALL(test_time_and_parameter_flow())
    ARIADNE_TEST_CALL(test_flow_derivatives())
    return ARIADNE_TEST_FAILURES;
}
