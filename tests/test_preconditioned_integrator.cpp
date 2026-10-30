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

#include <cmath>
#include <memory>
#include <sstream>
#include "integration/integrator.hpp"
#include "function/taylor_function.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
Void check_factorisation(PreconditionedTaylorSeriesState const& state,
                         ValidatedVectorMultivariateFunctionPatch const& original,
                         Vector<FloatDPBounds> const& point) {
    auto expected=original(point);
    auto local=state.normalised_mapping()(point);
    ARIADNE_TEST_EQUALS(state.parameter_domain().dimension(),original.argument_size())
    ARIADNE_TEST_EQUALS(state.local_domain().dimension(),original.result_size())
    for(SizeType i=0u;i<expected.size();++i) {
        FloatDPBounds reconstructed(state.centre()[i]);
        for(SizeType j=0u;j<local.size();++j) {
            reconstructed+=FloatDPBounds(state.linear_map()[i][j])*local[j];
        }
        ARIADNE_TEST_ASSERT(refines(expected[i],reconstructed))
    }
}

Void test_preconditioning_geometry() {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain={{-1.0_x,1.0_x},{-1.0_x,1.0_x}};
    auto const& factory=integrator.function_factory();
    auto id=factory.create_identity(domain);
    auto mapping=factory.create_zeros(2u,domain);
    mapping[0u]=FloatDPBounds(1,dp)+id[0u]+FloatDPBounds(3,dp)*id[1u];
    mapping[1u]=FloatDPBounds(2,dp)+FloatDPBounds(2,dp)*id[0u]+hlf(id[1u]);
    for(auto mode:{TaylorSeriesPreconditioning::IDENTITY,TaylorSeriesPreconditioning::QR}) {
        integrator.set_preconditioning(mode);
        auto state=integrator.precondition(mapping);
        check_factorisation(state,mapping,Vector<FloatDPBounds>({1,-1},dp));
        check_factorisation(state,mapping,Vector<FloatDPBounds>({0,0},dp));
        for(SizeType i=0u;i<2u;++i) {
            for(SizeType j=0u;j<2u;++j) {
                double dot=0.0;
                for(SizeType k=0u;k<2u;++k) {
                    dot+=state.linear_map()[k][i].get_d()*state.linear_map()[k][j].get_d();
                }
                ARIADNE_TEST_WITHIN(dot,(i==j ? 1.0 : 0.0),1e-12)
            }
        }
        std::ostringstream text;
        text << static_cast<IntegratorInterface const&>(integrator);
        ARIADNE_TEST_ASSERT(text.str().find(mode==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY")!=std::string::npos)
        std::unique_ptr<PreconditionedGradedTaylorSeriesIntegrator> copy(integrator.clone());
        ARIADNE_TEST_ASSERT(copy->preconditioning()==mode)
    }
    // Extra enclosure parameters intentionally retain identity orientation.
    ExactBoxType parameter_domain(3u,ExactIntervalType(-1.0_x,1.0_x));
    auto parameters=factory.create_identity(parameter_domain);
    auto rectangular=factory.create_zeros(2u,parameter_domain);
    rectangular[0u]=parameters[0u]+parameters[2u];
    rectangular[1u]=parameters[1u];
    check_factorisation(integrator.precondition(rectangular),rectangular,
                        Vector<FloatDPBounds>({1,-1,1},dp));
}

Void test_rank_deficient_geometry() {
    PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
    ExactBoxType domain(2u,ExactIntervalType(-1.0_x,1.0_x));
    auto const& factory=integrator.function_factory();
    auto constant=factory.create_zeros(2u,domain);
    constant[0u]=factory.create_constant(domain,1);
    constant[1u]=factory.create_constant(domain,2);
    auto state=integrator.precondition(constant);
    check_factorisation(state,constant,Vector<FloatDPBounds>({0,0},dp));
    ARIADNE_TEST_EQUALS(state.linear_map()[0u][0u],FloatDP(1,dp))
    ARIADNE_TEST_EQUALS(state.linear_map()[1u][1u],FloatDP(1,dp))
}

Void check_rotation_endpoint(PreconditionedTaylorSeriesStep const& step,StepSizeType const& elapsed) {
    Vector<FloatDPBounds> point({1,1},dp);
    auto value=step.evolved_mapping()(point);
    FloatDPBounds t(elapsed,dp);
    Vector<FloatDPBounds> expected({cos(t)-sin(t),sin(t)+cos(t)});
    for(SizeType i=0u;i<2u;++i) {
        ARIADNE_TEST_ASSERT(refines(expected[i],value[i]))
    }
    ARIADNE_TEST_EQUALS(step.final_state().centre().size(),2u)
    ARIADNE_TEST_EQUALS(step.final_state().parameter_domain().dimension(),2u)
    ARIADNE_TEST_ASSERT(definitely(is_bounded(step.flowpipe_mapping().range())))
}

Void test_rotating_flow_and_carried_state() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(2u,0u);
    auto y=EffectiveScalarMultivariateFunction::coordinate(2u,1u);
    ValidatedVectorMultivariateFunction field=EffectiveVectorMultivariateFunction({-y,x});
    ExactBoxType domain(2u,ExactIntervalType(0.875_x,1.125_x));
    StepSizeType h(0.03125_x);
    for(auto mode:{TaylorSeriesPreconditioning::IDENTITY,TaylorSeriesPreconditioning::QR}) {
      for(Bool diagnostics:{false,true}) {
        ARIADNE_PRINT_TEST_CASE_TITLE("Two rotation steps, mode " << (mode==TaylorSeriesPreconditioning::QR ? "QR" : "IDENTITY") << ", diagnostics " << diagnostics)
        PreconditionedGradedTaylorSeriesIntegrator integrator(step_maximum_error=1e-6);
        integrator.set_preconditioning(mode);
        integrator.set_diagnostics(diagnostics);
        integrator.set_carried_expansion_diagnostics(diagnostics);
        ARIADNE_TEST_EQUALS(integrator.diagnostics(),diagnostics)
        ARIADNE_TEST_EQUALS(integrator.carried_expansion_diagnostics(),diagnostics)
        auto initial=integrator.function_factory().create_identity(domain);
        auto first=integrator.step(field,integrator.precondition(initial),suggest(h));
        ARIADNE_TEST_COMPARE(first.time_step(),>,StepSizeType(0))
        ARIADNE_TEST_COMPARE(first.time_step(),<=,h)
        check_rotation_endpoint(first,first.time_step());
        auto second=integrator.step(field,first.final_state(),suggest(h));
        check_rotation_endpoint(second,first.time_step()+second.time_step());
      }
    }
}
}

int main() {
    ARIADNE_TEST_CALL(test_preconditioning_geometry())
    ARIADNE_TEST_CALL(test_rank_deficient_geometry())
    ARIADNE_TEST_CALL(test_rotating_flow_and_carried_state())
    return ARIADNE_TEST_FAILURES;
}
