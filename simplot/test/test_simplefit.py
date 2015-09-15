
import itertools
import math
import random
import string
import unittest

import numpy as np

from simplot.mc.montecarlo import ToyMC
from simplot.mc.generators import GaussianGenerator
from simplot.mc.statistics import Mean, StandardDeviation, calculate_statistics_from_toymc
from simplot.mc.likelihood import EventRateLikelihood, SumLikelihood
from simplot.mc.priors import GaussianPrior, CombinedPrior, OscillationParametersPrior
from simplot.binnedmodel.sample import Sample, BinnedSample, BinnedSampleWithOscillation, CombinedBinnedSample
from simplot.binnedmodel.systematics import Systematics, SplineSystematics, FluxSystematics, FluxAndSplineSystematics

from simplot.binnedmodel.simplemodel import SimpleMcBuilder

################################################################################

class TestSimpleFit(unittest.TestCase):

    def _buildtestmc(self, cachestr=None):
        systematics = [("x", [-5.0, 0.0, 5.0]),
                       ("y", [-5.0, 0.0, 5.0])]
        systematics = SplineSystematics(systematics)
        def gen(N):
            for _ in xrange(N):
                coord = np.random.poisson(size=2)
                yield coord, 1.0, [(-4.0, 1.0, 5.0), (-4.0, 1.0, 5.0)]
        binning = [("a", np.arange(0.0, 5.0)), ("b", np.arange(0.0, 5.0))]
        observables = ["a"]
        model = BinnedSample("simplemodel", binning, observables, gen(10**4), systematics=systematics, cache_name=cachestr)
        generator = GaussianGenerator(["x", "y"], [1.0, 2.0], [1.0, 2.0])
        toymc = ToyMC(model, generator)
        return toymc

    def test_build_simple_model(self):
        toymc1 = self._buildtestmc()
        toymc2 = SimpleMcBuilder().build("testmodel", toymc1)
        return

    def test_eval_model(self):
        toymc1 = self._buildtestmc()
        toymc2 = SimpleMcBuilder().build("testmodel", toymc1)
        stat1 = [Mean(), StandardDeviation()]
        stat2 = [Mean(), StandardDeviation()]
        npe = 10**4
        calculate_statistics_from_toymc(toymc1, stat1, npe=npe)
        calculate_statistics_from_toymc(toymc2, stat2, npe=npe)
        for st1, st2 in zip(stat1, stat2):
            v1, e1 = st1.eval(), st1.err()
            v2, e2 = st2.eval(), st2.err()
            for v1, e1, v2, e2 in zip(np.nditer(v1), np.nditer(e1), np.nditer(v2), np.nditer(e2)):
                delta = 5.0 * np.sqrt(e1**2 + e2**2)
                self.assertAlmostEquals(v1, v2, delta=delta)
        return

################################################################################

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()

