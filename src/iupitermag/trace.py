import numpy as np

import iupitermag._core as _iu

from . import currentsheet, internal


def trace_field_to_planet(
    start_positions: np.ndarray,
    internal_field: str | internal.InternalField = "JRM33",
    currentsheet_field: str | currentsheet.CurrentSheetField = "CON2020",
):
    """
    Trace the magnetic field from a collection of points to the planet.

    Args:
        start_positions (np.ndarray): Array of shape (N, 3) where N is the number
            of points and indices [:,0], [:, 1], and [:, 2] represent the X, Y,
            and Z coordinates of each point in the IAU coordinate system.
        internal_field (str | internal.InternalField): The internal field to use (default="JRM33").
        currentsheet_field (str | currentsheet.CurrentSheetField): The current sheet field to use
            (default="CON2020").

    Returns:
        traces (list[np.ndarray]): List of traces where each trace is of shape (M, 3) where M is
            the number of coordinates of each trace.  Each point of the trace is described by its
            cartesian representation in the IAU coordinate system.  In total there should be N
            traces, one for each starting position.
    """
    if isinstance(internal_field, str):
        internal_field = internal.InternalField(internal_field)

    if isinstance(currentsheet_field, str):
        currentsheet_field = currentsheet.CurrentSheetField(currentsheet_field)

    return _iu.trace_field_to_planet(
        np.asarray(start_positions, dtype=float).reshape(-1, 3),
        internal_field._field,
        currentsheet_field._field,
    )
