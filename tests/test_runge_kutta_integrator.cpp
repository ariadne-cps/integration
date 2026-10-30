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

#include "integration/runge_kutta_integrator.hpp"
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {
Void test_classical_fourth_order_step() {
    auto x=EffectiveScalarMultivariateFunction::coordinate(1u,0u);
    ApproximateVectorMultivariateFunction f=EffectiveVectorMultivariateFunction({x});
    FloatDPApproximationVector initial({FloatDPApproximation(1,dp)});
    RungeKutta4Integrator integrator(0.5);
    auto next=integrator.step(f,initial,FloatDPApproximation(0.5,dp));
    // For x'=x the fourth-order Taylor polynomial at h=1/2 is 211/128.
    ARIADNE_TEST_WITHIN(next[0u].get_d(),211.0/128.0,1e-14)
    ARIADNE_TEST_EQUALS(initial[0u].get_d(),1.0)
}

Void test_evolution_uses_each_instance_step() {
    auto one=EffectiveScalarMultivariateFunction::constant(1u,1);
    ApproximateVectorMultivariateFunction f=EffectiveVectorMultivariateFunction({one});
    FloatDPApproximationVector initial({FloatDPApproximation(2,dp)});
    RungeKutta4Integrator coarse(0.25);
    RungeKutta4Integrator fine(0.125);
    auto coarse_path=coarse.evolve(f,initial,FloatDPApproximation(0.5,dp));
    auto fine_path=fine.evolve(f,initial,FloatDPApproximation(0.5,dp));
    ARIADNE_TEST_EQUALS(coarse_path.size(),3u)
    ARIADNE_TEST_EQUALS(fine_path.size(),5u)
    ARIADNE_TEST_EQUALS(coarse_path[1u].first.get_d(),0.25)
    ARIADNE_TEST_EQUALS(fine_path[1u].first.get_d(),0.125)
    ARIADNE_TEST_EQUALS(coarse_path.back().second[0u].get_d(),2.5)
    ARIADNE_TEST_EQUALS(fine_path.back().second[0u].get_d(),2.5)
    auto zero_path=fine.evolve(f,initial,FloatDPApproximation(0,dp));
    ARIADNE_TEST_EQUALS(zero_path.size(),1u)
    ARIADNE_TEST_EQUALS(zero_path[0u].second[0u].get_d(),2.0)
}
}

int main() {
    ARIADNE_TEST_CALL(test_classical_fourth_order_step())
    ARIADNE_TEST_CALL(test_evolution_uses_each_instance_step())
    return ARIADNE_TEST_FAILURES;
}
