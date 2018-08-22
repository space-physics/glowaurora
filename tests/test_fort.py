from glowaurora import glowfort
from pytest import approx
import pytest


def test_glow():
    # electron precipitation
    # First enact "glow" subroutine, which calls QBACK, ETRANS and GCHEM among others

    glowfort.glow()  # no args

    # %% ver and constituants
    """
    using common block CGLOW, instead use new GLOW for module
    """
    zceta = glowfort.cglow.zceta.T
    zeta = glowfort.cglow.zeta.T[:, :11]
    zcsum = zceta.sum(axis=-1)[:, :11]
    assert zcsum == approx(zeta)


if __name__ == '__main__':
    pytest.main(['-xv', __file__])
