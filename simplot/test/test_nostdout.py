from simplot.nostdout import nostdout

import unittest
import sys
from StringIO import StringIO

_MSG_ERR = "TEST stderr message"
_MSG_OUT = "TEST stdout message"

def _test_stdout():
    print _MSG_OUT,

def _test_stderr():
    print >>sys.stderr, _MSG_ERR,

def _runtest(*args, **kwargs):
        strio = StringIO()
        stdout = sys.stdout
        stderr = sys.stderr
        sys.stdout = strio
        sys.stderr = strio
        with nostdout(*args, **kwargs):
            _test_stdout()
            _test_stderr()
        sys.stdout = stdout
        sys.stderr = stderr
        return strio.getvalue()

class TestNoStdOut(unittest.TestCase):
    def test_nostdout(self):
        self.assertEquals(_runtest(stdout=True, stderr=False), _MSG_ERR)
    def test_nostderr(self):
        self.assertEquals(_runtest(stdout=False, stderr=True), _MSG_OUT)
    def test_nostdboth(self):
        self.assertEquals(_runtest(stdout=True, stderr=True), "")

def _debugstdout():
    print "Messages ON"
    _test_stdout()
    _test_stderr()
    print "Messages OFF"
    with nostdout(stdout=True, stderr=False):
        _test_stdout()
    with nostdout(stdout=False, stderr=True):
        _test_stderr()
    print "Messages ON"
    _test_stdout()
    _test_stderr()
    print "done."
    return

def main():
    unittest.main()

if __name__ == "__main__":
    main()
