
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
