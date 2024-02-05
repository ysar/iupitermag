import numpy as np
import iupitermag


def test_planet_field():

    x = 10.
    y = 0.
    z = 0.

    inArr = np.array([[x, y, z]])

    outArr = iupitermag.get_planet_field(inArr, "JRM09", True, True, 10)
    
    expectedArr = np.array([[-131.3754, -24.2054 , -400.6375]])

    assert(np.allclose(outArr, expectedArr, rtol=1e-3))