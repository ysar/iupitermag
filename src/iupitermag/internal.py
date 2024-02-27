import iupitermag.iupitermag as _iu

class InternalField:

    def __init__(self, typefield="", g=None, h=None, degree=None):
        """
        Class for handling the planet's internal field. 

        Use: Initialize using 
            InternalField()

        Args:
            typefield (str): Type of planet field. Allowed values are 
                "JRM09", "JRM33", "Custom"
            g (np.array): Legendre coefficient array g[n ,m]
            h (np.array): Legendre coefficient array h[n, m]
            degree (int): Degree of field

            g, h, and degree do not need to be specified for non-"Custom" type.

        Returns:
            InternalField class object
        """
        self._field = _iu.PyInternalField(typefield, g, h, degree)


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


    def get_coefficients(self):
        """
        Get the coeffiicients of the defined internal field. 

        Returns:
            (g, h) (np.array, np.array): Legendre coefficients in units of nT
        """
        return self._field.get_coefficients()
    
    
    def map_calc_field(self, positions):
        return self._field.map_calc_field(positions)
    
    def parmap_calc_field(self, positions):
        return self._field.parmap_calc_field(positions)
