import cProfile
import pstats
import argparse

###############################################################################

def run_profiler(func, *args, **kwargs):
    profile = cProfile.Profile()
    profile.enable()
    func(*args, **kwargs)
    profile.disable()
    return profile

###############################################################################

def profile_func(func, args=[], kwargs={}, num=20):
    prof = run_profiler(func, *args, **kwargs)
    ps = pstats.Stats(prof)
    for key in ["time", "cumulative"]:
        print "--- top {} sorted by {}".format(num, key)
        ps.sort_stats(key).print_stats(num)
    return

###############################################################################

def main():
    def sum_inv(N=10000, start=1):
        ret = 0.0
        print start
        for x in xrange(start, N):
            ret += 1.0/float(x)
        return ret
    profile_func(sum_inv)
    profile_func(sum_inv, [100])
    profile_func(sum_inv, [], {"N":100, "start":10})
    profile_func(sum_inv, [100], {"start":10})
    return

###############################################################################

if __name__ == "__main__":
    main()
