import numpy as np

import iupitermag


def test_planet_field_point():

    x = 10.0
    y = 0.0
    z = 0.0

    in_arr = np.array([np.sqrt(x**2 + y**2 + z**2), np.arccos(z / x), np.arctan2(y, x)])

    internal_field = iupitermag.InternalField("JRM09", degree=10)

    out_arr = internal_field.calc_field(in_arr[0], in_arr[1], in_arr[2])
    expectedArr = np.array([-131.3754, 400.6375, -24.2054])

    assert np.allclose(out_arr, expectedArr, rtol=1e-3)
