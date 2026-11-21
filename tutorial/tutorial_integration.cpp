/***************************************************************************
 *            tutorial_integration.cpp
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

/** @file tutorial/tutorial_integration.cpp
 *  @brief Demonstrate the Integration module's standalone C++ API.
 */

#include <iostream>
#include <string>

#include "ariadne-integration.hpp"

using namespace Ariadne;
using std::cout;

namespace {

void section(const std::string& title) {
    cout << "\n== " << title << " ==\n";
}

template<class T>
void show(const std::string& label, const T& value) {
    cout << label << ": " << value << "\n";
}

} // namespace

int main() {
    section("Validated ODE integration with Taylor-Picard");

    TaylorPicardIntegrator integrator(1e-8);
    integrator.set_minimum_temporal_order(2);
    integrator.set_maximum_temporal_order(16);
    integrator.set_step_maximum_error(1e-8);

    show("minimum temporal order", integrator.minimum_temporal_order());
    show("maximum temporal order", integrator.maximum_temporal_order());
    show("step maximum error", integrator.step_maximum_error());

    auto state = EffectiveScalarMultivariateFunction::coordinate(1, 0);
    EffectiveVectorMultivariateFunction logistic({
        state * (1 - state)
    });
    ExactBoxType initial_states({
        {0.25_x, 0.5_x}
    });
    StepSizeType step = 0.125_dy;

    FlowStepModelType flow =
        integrator.flow_step(logistic, initial_states, step);

    show("vector field", logistic);
    show("initial-state box", initial_states);
    show("requested step", step);
    show("validated flow patch", flow);
    show("flow domain", flow.domain());
    show("flow range", flow.range());

    section("Common validated step interface and continuation");
    const ExactBoxType initial({{0.0_x, 0.125_x}});
    const UpperBoxType bound({{-1.0_x, 1.0_x}});
    const StepSizeType h1(0.03125_x), h2(0.015625_x);
    const FlowStepIntegratorAdapter flow_adapter(integrator);
    const ValidatedStepIntegratorInterface& equation_solver = flow_adapter;
    const EffectiveVectorMultivariateFunction constant = {
        EffectiveScalarMultivariateFunction::constant(1u, 1)};
    const auto first = equation_solver.step(ValidatedIntegrationRequest(constant, initial, h1, bound));
    const auto second = equation_solver.step(ValidatedIntegrationRequest(constant, first, h2, bound));
    show("absolute final time", second.final_time());
    show("validated endpoint", second.evolve());
    show("structured endpoint domain", second.evolve_model().domain());
    show("reach sections at the final time", second.reach_sections(second.final_time()).size());

    section("Differential inclusion with one configured reach/evolve pair");
    const auto u = EffectiveFormula::coordinate(1u);
    const EffectiveVectorMultivariateFunction input_field =
        EffectiveVectorMultivariateFormulaFunction(2u, Vector<EffectiveFormula>{u});
    // The numerical field and solver are named objects. Requests borrow the
    // former and the adapter borrows the latter; both must remain alive.
    ValidatedInputFieldBoxAdapterDP inclusion(input_field,
        BoxDomainType({{-1.0_x, 1.0_x}}), BoxDomainType({{-1.0_x, 1.0_x}}));
    SelectedParametricMultifunctionIntegrator selected(integrator,
        InputApproximation(ZeroApproximation()), InputApproximation(ConstantApproximation()));
    const StructuredInclusionStepIntegratorAdapter inclusion_adapter(selected);
    const ValidatedStepIntegratorInterface& inclusion_solver = inclusion_adapter;
    const auto inclusion_first = inclusion_solver.step(ValidatedIntegrationRequest(inclusion, initial, h1, bound));
    const auto inclusion_second = inclusion_solver.step(ValidatedIntegrationRequest(inclusion, inclusion_first, h2, bound));
    show("inclusion endpoint", inclusion_second.evolve());
    show("retained state and input domain", inclusion_second.evolve_model().domain());

    section("An unsupported capability fails explicitly");
    try {
        equation_solver.step(ValidatedIntegrationRequest(inclusion, initial, h1, bound));
        return 1; // Unexpected: an equation solver cannot consume this field.
    } catch (const std::invalid_argument& error) {
        show("rejected request", error.what());
    }

    return 0;
}
