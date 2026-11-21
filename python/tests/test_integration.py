import gc
import unittest

from pyariadne_integration import *


def test_integration():
    integrator = TaylorPicardIntegrator(1e-8)
    integrator.minimum_temporal_order()
    integrator.maximum_temporal_order()


# unittest also permits running these regressions without installing pytest;
# the repository's pytest runner discovers the same cases in CI.
def exact_box(bounds):
    return FloatDPExactBox([{exact(lower): exact(upper)} for lower, upper in bounds])


class TestValidatedSteps(unittest.TestCase):
    def assert_encloses(self, box, lower, upper):
        expected = FloatDPBounds(exact(lower), exact(upper), dp)
        actual = FloatDPBounds(box[0].lower_bound(), box[0].upper_bound())
        self.assertTrue(refines(expected, actual))

    def test_adaptive_request_preserves_origin_and_mapping(self):
        field = EffectiveVectorMultivariateFunction([
            EffectiveScalarMultivariateFunction.constant(1, 1)])
        initial = exact_box([(0.0, 0.125)])
        maximum, origin = Dyadic("0.03125"), Dyadic(2)
        request = ValidatedIntegrationRequest.with_maximum_step(field, initial, maximum, origin)
        self.assertTrue(request.step_size_is_maximum())
        self.assertFalse(request.has_bounding_box())
        with self.assertRaises(RuntimeError):
            request.bounding_box()
        solver = FlowStepIntegratorAdapter(TaylorPicardIntegrator(1e-6))
        first = solver.step(request)
        self.assertEqual(first.initial_time(), origin)
        self.assertLessEqual(first.time_step(), maximum)
        second = solver.step(ValidatedIntegrationRequest.with_maximum_step(field, first, maximum))
        self.assertEqual(second.initial_time(), first.final_time())
        self.assertEqual(second.evolve_model().domain(), BoxDomainType(initial))
        self.assert_encloses(second.evolve(), 0.0625, 0.1875)

    def test_equation_request_survives_temporary_function_and_solver(self):
        initial = exact_box([(0.0, 0.125)])
        bound = FloatDPUpperBox(exact_box([(-1.0, 1.0)]))
        h = Dyadic("0.03125")
        request = ValidatedIntegrationRequest(
            EffectiveVectorMultivariateFunction([
                EffectiveScalarMultivariateFunction.constant(1, 1)]),
            initial, h, bound)
        solver = FlowStepIntegratorAdapter(TaylorSeriesIntegrator(1e-6, 8))
        gc.collect()
        first = solver.step(request)
        self.assert_encloses(first.evolve(), 0.03125, 0.15625)
        field = EffectiveVectorMultivariateFunction([
            EffectiveScalarMultivariateFunction.constant(1, 1)])
        second = solver.step(ValidatedIntegrationRequest(field, first, h, bound))
        self.assertEqual(second.final_time(), Dyadic("0.0625"))
        self.assert_encloses(second.evolve(), 0.0625, 0.1875)
        self.assertEqual(len(second.reach_sections(second.final_time())), 1)
        model = second.evolve_model()
        del first, second, solver
        gc.collect()
        self.assertEqual(model.domain(), BoxDomainType(initial))

    def test_structured_inclusion_retains_fields_solvers_and_input_choices(self):
        u = EffectiveScalarMultivariateFunction.coordinate(2, 1)
        initial = exact_box([(0.0, 0.125)])
        bound = FloatDPUpperBox(exact_box([(-1.0, 1.0)]))
        h1, h2 = Dyadic("0.03125"), Dyadic("0.015625")
        # Both the field and the configured backend are Python temporaries.
        # Their borrowing wrappers must retain them until no longer needed.
        request = ValidatedIntegrationRequest(
            ValidatedInputFieldBoxAdapterDP(
                EffectiveVectorMultivariateFunction([u]),
                BoxDomainType(exact_box([(-1.0, 1.0)])), BoxDomainType(exact_box([(-1.0, 1.0)]))),
            initial, h1, bound)
        solver = StructuredInclusionStepIntegratorAdapter(
            SelectedParametricMultifunctionIntegrator(
                TaylorSeriesIntegrator(1e-6, 8),
                InputApproximation(ZeroApproximation()),
                InputApproximation(ConstantApproximation())))
        gc.collect()
        first = solver.step(request)
        field = ValidatedInputFieldBoxAdapterDP(
            EffectiveVectorMultivariateFunction([u]),
            BoxDomainType(exact_box([(-1.0, 1.0)])), BoxDomainType(exact_box([(-1.0, 1.0)])))
        second = solver.step(ValidatedIntegrationRequest(field, first, h2, bound))
        self.assertEqual(second.final_time(), Dyadic("0.046875"))
        self.assert_encloses(second.evolve(), -0.046875, 0.171875)
        self.assertGreater(second.evolve_model().argument_size(), first.evolve_model().argument_size())
        self.assertTrue(second.has_reach_models())
        self.assertGreaterEqual(len(second.reach_sections(second.final_time())), 1)

    def test_timed_inclusion_parameters_and_absolute_time(self):
        t = EffectiveScalarMultivariateFunction.coordinate(4, 1)
        p = EffectiveScalarMultivariateFunction.coordinate(4, 2)
        u = EffectiveScalarMultivariateFunction.coordinate(4, 3)
        field = ValidatedTimedInputFieldBoxAdapterDP(
            EffectiveVectorMultivariateFunction([t + p + u]),
            BoxDomainType(exact_box([(-1.0, 1.0)])),
            IntervalDomainType(exact(0.0), exact(1.0)),
            BoxDomainType(exact_box([(0.125, 0.25)])),
            BoxDomainType(exact_box([(-0.125, 0.125)])))
        initial = exact_box([(0.0, 0.0)])
        bound = FloatDPUpperBox(exact_box([(-1.0, 1.0)]))
        parameters = exact_box([(0.125, 0.25)])
        h = Dyadic("0.03125")
        box_solver = MultiflowStepIntegratorAdapter(FirstOrderMultifunctionIntegrator())
        structured_solver = StructuredInclusionStepIntegratorAdapter(
            SelectedParametricMultifunctionIntegrator(TaylorSeriesIntegrator(1e-6, 8),
                InputApproximation(ZeroApproximation()), InputApproximation(ConstantApproximation())))
        for solver in (box_solver, structured_solver):
            first = solver.step(ValidatedIntegrationRequest(
                field, initial, Dyadic("0.5"), h, parameters, bound))
            second = solver.step(ValidatedIntegrationRequest(field, first, h, bound))
            self.assertEqual(second.initial_time(), Dyadic("0.53125"))
            self.assertEqual(second.final_time(), Dyadic("0.5625"))
            self.assertTrue(second.has_explicit_time())
            self.assertEqual(second.constant_parameter_domain(), parameters)
            # Exact integral of t+p+u(t) for t in [1/2,9/16].
            self.assert_encloses(second.evolve(), 0.033203125, 0.056640625)
        self.assertTrue(second.has_evolve_model())

    def test_unsupported_capability_and_failed_certification_raise(self):
        u = EffectiveScalarMultivariateFunction.coordinate(2, 1)
        field = ValidatedInputFieldBoxAdapterDP(
            EffectiveVectorMultivariateFunction([u]),
            BoxDomainType(exact_box([(-1.0, 1.0)])), BoxDomainType(exact_box([(-1.0, 1.0)])))
        initial = exact_box([(0.0, 0.125)])
        bound = FloatDPUpperBox(exact_box([(-1.0, 1.0)]))
        request = ValidatedIntegrationRequest(field, initial, Dyadic("0.03125"), bound)
        solver = FlowStepIntegratorAdapter(TaylorSeriesIntegrator(1e-6, 8))
        with self.assertRaises(ValueError):
            solver.step(request)
        with self.assertRaises(ValueError):
            SelectedParametricMultifunctionIntegrator(TaylorSeriesIntegrator(1e-6, 8),
                InputApproximation(AffineApproximation()), InputApproximation(ConstantApproximation()))
        with self.assertRaises(RuntimeError):
            MultiflowStepIntegratorAdapter(FirstOrderMultifunctionIntegrator()).step(
                ValidatedIntegrationRequest(field, initial, Dyadic("0.03125"), FloatDPUpperBox(initial)))
