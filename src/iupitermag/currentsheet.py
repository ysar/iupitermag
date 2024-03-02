import iupitermag.iupitermag as _iu

from .field import Field

class CurrentSheetField(Field):

    def __init__(self, typefield, params={}):
        self._field = _iu.PyCurrentSheetField(typefield, params)