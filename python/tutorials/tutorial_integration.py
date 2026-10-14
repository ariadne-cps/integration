#!/usr/bin/python3

from pyariadne import *


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


if __name__ == "__main__":
    tutorial_integration()
