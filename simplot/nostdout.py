import contextlib
import sys

@contextlib.contextmanager
def nostdout(stdout=True, stderr=True):
    """Suppress stdout and stderr.

    Example:
    >>> with nostdout():
    >>>     print "should not print"
    >>>

    stolen from: http://stackoverflow.com/questions/1809958/hide-stderr-output-in-unit-tests
    """
    savestdout = sys.stdout
    savestderr = sys.stderr
    class Devnull(object):
        def write(self, _): pass
        def flush(self): pass
    if stdout:
        sys.stdout = Devnull()
    if stderr:
        sys.stderr = Devnull()
    try:
        yield
    finally:
        sys.stdout = savestdout
        sys.stderr = savestderr
