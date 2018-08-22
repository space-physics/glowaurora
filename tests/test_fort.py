#!/usr/bin/env python
import glowaurora as ga
from pytest import approx
import pytest


def test_glow():
    # electron precipitation
    # First enact "glow" subroutine, which calls QBACK, ETRANS and GCHEM among others

    ga.glowfort.glow()  # no args

    # %% ver and constituants
    """
    using common block CGLOW, instead use new GLOW for module
    """
    zceta = ga.glowfort.cglow.zceta.T
    zeta = ga.glowfort.cglow.zeta.T[:, :11]
    zcsum = zceta.sum(axis=-1)[:, :11]
    assert zcsum == approx(zeta)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
