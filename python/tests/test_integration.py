from pyariadne import *


def test_integration():
    integrator = TaylorPicardIntegrator(1e-8)
    integrator.minimum_temporal_order()
    integrator.maximum_temporal_order()
