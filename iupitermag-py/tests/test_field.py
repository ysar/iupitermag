import numpy as np

import iupitermag


def test_internal_field():

    x = 10.0
    y = 0.0
    z = 0.0

    pos = np.array([np.sqrt(x**2 + y**2 + z**2), np.arccos(z / x), np.arctan2(y, x)])

    internal_field = iupitermag.InternalField("JRM09", degree=10)

    b_calc = internal_field.calc_field(pos[0], pos[1], pos[2])
    b_expected = np.array([-131.3754, 400.6375, -24.2054])
    assert np.allclose(b_expected, b_calc, rtol=1e-3)


def test_currentsheet_field():
    x = 20.2356
    y = 1.31
    z = -6.51

    currentsheet_field = iupitermag.CurrentSheetField("CON2020")

    b_calc = currentsheet_field.calc_field_xyz(x, y, z)
    b_expected = np.array([-31.20, 1.53, 23.60])
    assert np.allclose(b_expected, b_calc, rtol=1e-3)
