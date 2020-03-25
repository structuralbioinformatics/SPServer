class Singleton(type):
    '''
    This \href{}{Singleton} design has been obtained from an
    \href{http://stackoverflow.com/q/6760685/2806632}{entry in StackOverflow}.
    '''
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

    @staticmethod
    def _list():
        for cls in Singleton._instances:
            yield Singleton._instances[cls]
