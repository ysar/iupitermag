from .iupitermag import PyInternalField

class InternalField:

    def __init__(self, typefield="", g=None, h=None, degree=None):
        """

        """
        self._field = PyInternalField(typefield, g, h, degree)

    def calc_field(self, r, theta, phi):
        return self._field.calc_internal_field(r, theta, phi)
    
    def get_coefficients(self):
        return self._field.get_coefficients()
    
    def set_coefficients(self, g, h):
        pass

    def map_calc_field(self, positions):
        pass

    def loop_calc_field(self, positions):
        pass