import numpy as np

import iupitermag


def test_trace():

    start_pos = np.array([[-10.0, 2.0, 3.0]])

    internal_field = iupitermag.InternalField("JRM33", degree=10)
    currentsheet_field = iupitermag.CurrentSheetField("CON2020")

    traces = iupitermag.trace_field_to_planet(
        start_pos, internal_field, currentsheet_field
    )

    # The point traced to the planet in the northen hemisphere.
    first_expected = np.array([-0.52819934, -0.01459706, 0.77033209])
    assert np.allclose(first_expected, traces[0][0, :], rtol=1e-3)
