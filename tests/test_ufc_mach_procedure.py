"""Mach-front loads follow UFC 3-340-02 §2-13.2 / Problems 2A-2 and 2A-10."""
import math
import warnings

import numpy as np
import pytest

from ufc_blast.core import blast_params as bp
from ufc_blast.core.blast_params import compute_facade, compute_mach_front

FT = 0.3048 / 0.45359237 ** (1 / 3)      # ft/lb^1/3 -> m/kg^1/3
PSI = 6.894757
IMP = PSI / 0.45359237 ** (1 / 3)        # psi-ms/lb^1/3 -> kPa-ms/kg^1/3


@pytest.fixture(autouse=True)
def _quiet():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bp.load_ufc_tables()
        yield


def test_example_2a2_ground_values():
    """UFC Example 2A-2: W=25,000 lb, RG=300 ft, Hc=90 ft (hand-read values)."""
    W = 25000 * 0.45359237
    RG, Hc = 300 * 0.3048, 90 * 0.3048
    W_cbrt = W ** (1 / 3)
    alpha = math.degrees(math.atan2(RG, Hc))
    assert alpha == pytest.approx(73.3, abs=0.05)
    T = bp._tables
    Pra = T["pra_ground"].lookup(x=alpha, family=Hc / W_cbrt)
    ira = T["ira_ground"].lookup(x=alpha, family=Hc / W_cbrt)
    assert Pra / PSI == pytest.approx(10.1, rel=0.10)
    assert ira / IMP == pytest.approx(9.2, rel=0.10)
    Zp, _ = bp._inverse_lookup(T["ps0"], Pra)
    Zi, _ = bp._inverse_lookup(T["is"], ira)
    assert Zp / FT == pytest.approx(7.8, rel=0.05)
    assert Zi / FT == pytest.approx(5.7, rel=0.06)
    t0 = T["t0"].lookup(Zi) * W_cbrt
    assert t0 == pytest.approx(45.3, rel=0.03)   # 1.55 ms/lb^1/3 * 25,000^1/3


def test_mach_region_uses_ground_values_not_triple_point():
    f = compute_facade(R=30.0, W=200.0, Hc=5.0, width=20.0, height=10.0, step=1.0)
    ref = compute_mach_front(R=30.0, dx=0.0, W=200.0, Hc=5.0)
    low = f.result_map[(0.0, -5.0)]
    assert low.Pr_alpha == pytest.approx(ref.Pr_alpha)
    assert low.ir_alpha == pytest.approx(ref.ir_alpha)


def test_mach_wave_stronger_than_free_air_incident():
    """UFC: Mach-front parameters are 'somewhat larger' than the incident wave."""
    W, R, Hc = 200.0, 30.0, 5.0
    m = compute_mach_front(R=R, dx=0.0, W=W, Hc=Hc)
    free = bp.compute_point(math.hypot(R, Hc), 0.0, W)
    assert m.Ps0 > free.Ps0


def test_normal_incidence_wall_reflection_uses_fig_2_193():
    m = compute_mach_front(R=30.0, dx=0.0, W=200.0, Hc=5.0)
    Ca0 = bp._tables["calpha"].lookup(x=0.0, family=m.Ps0)
    assert m.Pr_alpha == pytest.approx(Ca0 * m.Ps0)


def test_clipped_b_is_flagged():
    Pr, t0 = np.array([100.0]), np.array([0.02])
    assert bp._b_clipped(Pr, t0, np.array([0.99 * 100 * 0.02 / 2]))[0] == False  # noqa: E712
    assert bp._b_clipped(Pr, t0, np.array([1.2 * 100 * 0.02 / 2]))[0] == True    # noqa: E712


def test_threshold_warning():
    W = 100.0
    Hc = 0.4 * W ** (1 / 3)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        compute_facade(R=30.0, W=W, Hc=Hc, width=4.0, height=4.0, step=2.0)
    assert any("chart-selection rule" in str(x.message) for x in w)
