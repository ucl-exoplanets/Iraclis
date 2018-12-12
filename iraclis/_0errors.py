from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


class IraclisError(BaseException):
    pass


class IraclisLibraryError(IraclisError):
    pass


class IraclisFileError(IraclisError):
    pass


class IraclisProcessError(IraclisError):
    pass


class IraclisInputError(IraclisError):
    pass
