import iupitermag._core as _iu

from .field import Field


class CurrentSheetField(Field):
    def __init__(self, typefield, params={}, integration_type="analytic"):
        """
        Class for working with Jupiter current sheet models.

        Args:
            typefield (str): Type of current sheet field. Current available options are
                'CON2020' or 'Custom'. If using a named field, the `params` argument
                does not need to be used.

            params (dict): The parameters of the current sheet. Only necessary for 'Custom' field.

            integration_type (str): Type of integration used. Options are 'analytic' (default) or
                'integral' (not currently implemented).
        """
        self._field = _iu.PyCurrentSheetField(typefield, params, integration_type)

    def get_params(self):
        """
        Get the parameters of this current sheet field model.

        Returns:
            params (dict): Dictionary containing all parameters for this current sheet model.
        """
        return self._field.get_params()
