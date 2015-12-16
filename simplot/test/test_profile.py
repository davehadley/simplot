import unittest
from simplot.profile import profile_func
from simplot.nostdout import nostdout

def sum_inv(N=10000, start=1):
    ret = 0.0
    for x in xrange(start, N):
        ret += 1.0/float(x)
    return ret

class TestProfile(unittest.TestCase):
    def test_prof_noargs(self):
        with nostdout():
            profile_func(sum_inv)
    def test_prof_args(self):
        with nostdout():
            profile_func(sum_inv, [100])
    def test_prof_kwargs(self):
        with nostdout():
            profile_func(sum_inv, [], {"N":100, "start":10})
    def test_prof_both(self):
        with nostdout():
            profile_func(sum_inv, [100], {"start":10})

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
