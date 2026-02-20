import iupitermag._core as _iu

from .field import Field


class CurrentSheetField(Field):
    def __init__(self, typefield, params={}, integration_type="analytic"):
        self._field = _iu.PyCurrentSheetField(typefield, params, integration_type)
