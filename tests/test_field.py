import numpy as np
import iupitermag


def test_planet_field_point():

    x = 10.
    y = 0.
    z = 0.

    inArr = np.array([
        np.sqrt(x**2 + y**2 + z**2), 
        np.arccos(z / x),
        np.arctan2(y, x)
        ])

    internalField = iupitermag.InternalField("JRM09", degree=10)

    outArr = internalField.calc_field(inArr[0], inArr[1], inArr[2])
    expectedArr = np.array([-131.3754, 400.6375, -24.2054])
    
    assert(np.allclose(outArr, expectedArr, rtol=1e-3))


