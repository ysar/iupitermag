import numpy as np

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
            field (np.array): Magnetic field components [Br, Btheta, Bphi] in a
                array.
        """
        return self._field.calc_field(r, theta, phi)

    def calc_field_xyz(self, x, y, z):
        """
        Calculate the field at a location defined by (x, y, z) in a Cartesian
        coordinate system.

        Args:
            x (float): The X coordinate.
            y (float): The Y coordinate.
            z (float): The Z coordinate.

        Returns:
            field (np.ndarray): An array containing the magnetic field vector in
                cartesian coordinates (Bx, By, Bz)
        """
        return self._field.calc_field_xyz(x, y, z)

    def map_calc_field(self, positions):
        """
        Calculates [Br, Btheta, Bphi] for a collection of N points.

        Args:
            positions (np.ndarray): Array of spherical coordinates of shape (N, 3),
                where the last index refers to (r, theta, phi).

        Returns:
            field_rtp (np.ndarray): Array of (Br, Btheta, Bphi) of shape (N, 3)
                for each point.
        """
        return self._field.map_calc_field(
            np.asarray(positions, dtype=float).reshape(-1, 3)
        )

    def parmap_calc_field(self, positions):
        """
        Calculates [Br, Btheta, Bphi] for a collection of N points using Rayon
        for parallelization.

        Args:
            positions (np.ndarray): Array of spherical coordinates of shape (N, 3),
                where the last index refers to (r, theta, phi).

        Returns:
            field_rtp (np.ndarray): Array of (Br, Btheta, Bphi) of shape (N, 3)
                for each point.
        """
        return self._field.parmap_calc_field(
            np.asarray(positions, dtype=float).reshape(-1, 3)
        )

    def map_calc_field_xyz(self, positions):
        """
        Calculates [Bx, By, Bz] for a collection of N points. Uses cartesian coordinates
        as input.

        Args:
            positions (np.ndarray): Array of spherical coordinates of shape (N, 3),
                where the last index refers to (x, y, z) coordinates.

        Returns:
            field_xyz (np.ndarray): Array of shape (N, 3) containing (Bx, By, Bz)
                for each point.
        """
        return self._field.map_calc_field_xyz(
            np.asarray(positions, dtype=float).reshape(-1, 3)
        )

    def parmap_calc_field_xyz(self, positions):
        """
        Calculates [Bx, By, Bz] for a collection of N points. Uses cartesian coordinates
        as input and Rayon for parallelization.

        Args:
            positions (np.ndarray): Array of spherical coordinates of shape (N, 3),
                where the last index refers to (x, y, z) coordinates.

        Returns:
            field_xyz (np.ndarray): Array of shape (N, 3) containing (Bx, By, Bz)
                for each point.
        """
        return self._field.parmap_calc_field_xyz(
            np.asarray(positions, dtype=float).reshape(-1, 3)
        )
