import iupitermag._core as _iu

from .field import Field


class InternalField(Field):
    def __init__(self, typefield="", g=None, h=None, degree=None):
        """
        Class for handling the planet's internal field.

        Args:
            typefield (str): Type of planet field. Allowed values are
                "JRM09", "JRM33", "Custom"
            g (np.array): Legendre coefficient array g[n ,m]
            h (np.array): Legendre coefficient array h[n, m]
            degree (int): Degree of field

            g, h, and degree do not need to be specified for named fields.

        Returns:
            InternalField class object
        """
        self._field = _iu.PyInternalField(typefield, g, h, degree)

    def get_coefficients(self):
        """
        Get the coeffiicients of the defined internal field.

        Returns:
            g (np.array): Legendre coefficient (g) in units of nT.
            h (np.array): Legendre coefficient (h) in units of nT.
        """
        return self._field.get_coefficients()
