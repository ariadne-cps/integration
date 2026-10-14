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

    return 0;
}
