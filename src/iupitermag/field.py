import iupitermag._core as _iu


class Field:
    _field: _iu.PyInternalField | _iu.PyCurrentSheetField

    def __init__(self):
        """
        Base class for all field types. Contains generic routines like
        calc_field, map_calc_field, trace_field, etc.
        """
        pass

    def calc_field(self, r, theta, phi):
        """
        Calculate the field strength at a location defined in a spherical
        planetocentric coordinate system.

        Args:
            r (float): Radius in planetary radii
            theta (float): Co-latitude in radians
            phi (float): Azimuth in radians

        Returns:
            field (np.array): Magnetic field components [Br, Btheta, Bphi]
        """
        return self._field.calc_field(r, theta, phi)

    def calc_field_xyz(self, x, y, z):
        """
        Calculate the field strength at a location defined by (x, y, z) in a Cartesian
        coordinate system.
        """
        return self._field.calc_field_xyz(x, y, z)

    def get_params(self):
        """
        Get the parameters of the field. Varies for each field type.
        """
        return self._field.get_params()

    def map_calc_field(self, positions):
        return self._field.map_calc_field(positions)

    def parmap_calc_field(self, positions):
        return self._field.parmap_calc_field(positions)

    def map_calc_field_xyz(self, positions):
        return self._field.map_calc_field_xyz(positions)

    def parmap_calc_field_xyz(self, positions):
        return self._field.parmap_calc_field_xyz(positions)
