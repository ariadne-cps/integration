#!/usr/bin/python3

from pyariadne_integration import *


def section(title):
    print(f"\n== {title} ==")


def show(label, value):
    print(f"{label}: {value}")


def exact_box(bounds):
    return FloatDPExactBox([
        {exact(lower): exact(upper)}
        for lower, upper in bounds
    ])


def tutorial_integration():
    section("Validated ODE integration with Taylor-Picard")

    integrator = TaylorPicardIntegrator(1e-8)
    integrator.set_minimum_temporal_order(2)
    integrator.set_maximum_temporal_order(16)
    integrator.set_step_maximum_error(1e-8)

    show("integrator", integrator)
    show("minimum temporal order", integrator.minimum_temporal_order())
    show("maximum temporal order", integrator.maximum_temporal_order())
    show("step maximum error", integrator.step_maximum_error())

    state = EffectiveScalarMultivariateFunction.coordinate(1, 0)
    logistic = EffectiveVectorMultivariateFunction([
        state * (1 - state),
    ])
    initial_states = exact_box([
        (0.25, 0.5),
    ])
    step = Dyadic("0.125")

    flow = integrator.flow_step(logistic, initial_states, step)

    show("vector field", logistic)
    show("initial-state box", initial_states)
    show("requested step", step)
    show("validated flow patch", flow)
    show("flow domain", flow.domain())
    show("flow range", flow.range())

    section("Common validated step interface and two-step continuation")
    initial = exact_box([(0.0, 0.125)])
    bound = FloatDPUpperBox(exact_box([(-1.0, 1.0)]))
    h1, h2 = Dyadic("0.03125"), Dyadic("0.015625")
    solver = FlowStepIntegratorAdapter(integrator)
    constant = EffectiveVectorMultivariateFunction([
        EffectiveScalarMultivariateFunction.constant(1, 1)])
    first = solver.step(ValidatedIntegrationRequest(constant, initial, h1, bound))
    second = solver.step(ValidatedIntegrationRequest(constant, first, h2, bound))
    show("absolute final time", second.final_time())
    show("endpoint", second.evolve())
    show("endpoint domain", second.evolve_model().domain())
    show("reach sections", len(second.reach_sections(second.final_time())))

    section("Structured differential inclusion")
    u = EffectiveScalarMultivariateFunction.coordinate(2, 1)
    field = ValidatedInputFieldBoxAdapterDP(
        EffectiveVectorMultivariateFunction([u]),
        BoxDomainType(exact_box([(-1.0, 1.0)])),
        BoxDomainType(exact_box([(-1.0, 1.0)])))
    inclusion_solver = StructuredInclusionStepIntegratorAdapter(
        SelectedParametricMultifunctionIntegrator(
            integrator, InputApproximation(ZeroApproximation()),
            InputApproximation(ConstantApproximation())))
    first = inclusion_solver.step(ValidatedIntegrationRequest(field, initial, h1, bound))
    second = inclusion_solver.step(ValidatedIntegrationRequest(field, first, h2, bound))
    show("inclusion endpoint", second.evolve())
    show("retained state and input domain", second.evolve_model().domain())

    section("An unsupported capability fails explicitly")
    try:
        solver.step(ValidatedIntegrationRequest(field, initial, h1, bound))
    except ValueError as error:
        show("rejected request", error)
    else:
        raise AssertionError("An equation solver accepted a multifunction")


if __name__ == "__main__":
    tutorial_integration()
