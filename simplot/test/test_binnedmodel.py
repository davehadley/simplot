
import itertools
import math
import random
import string
import unittest

import numpy as np

from simplot.mc.montecarlo import ToyMC
from simplot.mc.statistics import Mean, StandardDeviation, calculate_statistics_from_toymc
from simplot.mc.likelihood import EventRateLikelihood, SumLikelihood
from simplot.mc.priors import GaussianPrior, CombinedPrior, OscillationParametersPrior
from simplot.binnedmodel.sample import BinnedSample, BinnedSampleWithOscillation, CombinedBinnedSample

class TestModel(unittest.TestCase):

    def test_model_building(self):
        #run twice to test caching
        cachestr = "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(25))
        for _ in xrange(2):
            self._buildmodelnoosc(cachestr=cachestr)

    def test_generate_mc(self):
        _, toymc, _ = self._buildmodelnoosc()
        npe = 100
        for _ in xrange(npe):
            toymc()
        return

    def test_model_building_withosc(self):
        #run twice to test caching
        cachestr = "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(25))
        for _ in xrange(2):
            self._buildmodelwithosc(cachestr=cachestr)

    def test_generate_mc_withosc(self):
        _, toymc, _ = self._buildmodelwithosc()
        npe = 100
        for _ in xrange(npe):
            toymc()
        return

    def _normalise(self, arr, norm=1.0):
        return arr * (norm/np.sum(arr))

    def test_generate_mc_values(self):
        _, toymc, _ = self._buildmodelnoosc()
        mean = Mean()
        stddev = StandardDeviation()
        statistics = [mean, stddev]
        calculate_statistics_from_toymc(toymc, statistics, 1000)
        #calculate the expected mean
        binedges = np.linspace(0.0, 5.0, num=100.0) 
        x = (binedges[1:] + binedges[:-1]) / 2.0
        mu = 1.0
        sigma = 0.25
        signal = self._normalise(1.0/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2) ), norm=0.5*10**5)
        background = self._normalise(np.ones(shape=signal.shape), norm=0.5*10**5)
        expectedmean = signal + background
        serr = 0.1
        berr = 0.1
        expectedsigma = np.sqrt(np.power(serr*signal, 2) + np.power(berr*background, 2)) * (mean.eval()/expectedmean) # correct sigma for stat err on mean
        #check prediction
        #import matplotlib.pyplot as plt
        #plt.errorbar(x, mean.eval(), yerr=np.sqrt(mean.eval()), color="blue")
        #plt.scatter(x, expectedmean, color="red")
        #plt.show()
        #raw_input("wait")
        for m, e, mex in itertools.izip_longest(mean.eval(), mean.err(), expectedmean):
            e = np.sqrt(e**2 + mex)
            self.assertAlmostEquals(m, mex, delta=5.0*e)
        for s, e, sex in itertools.izip_longest(stddev.eval(), stddev.err(), expectedsigma):
            self.assertAlmostEquals(s, sex, delta=3.0*e)
        return

    def _buildmodelnoosc(self, cachestr="testnoosc"):
        return self._buildmodel(withosc=False, cachestr=cachestr)

    def _buildmodelwithosc(self, cachestr="testwithosc"):
        return self._buildmodel(withosc=True, cachestr=cachestr)

    def _buildmodel(self, withosc, cachestr="testnoosc"):
        #build a simple dummy model
        binning = [("reco_energy", np.linspace(0.0, 5.0, num=100.0)),
                       ("true_energy", np.linspace(0.0, 5.0, num=25.0)),
                        ("true_nupdg", [0.0, 1.0, 2.0, 3.0, 4.0]),
        ]
        signalsyst = np.array([[-4.0, 1.0, 6.0], [1.0, 1.0, 1.0]])
        bkgdsyst = np.array([[1.0, 1.0, 1.0], [-4.0, 1.0, 6.0]])
        def gen(N):
            for _ in xrange(N):
                nupdg = np.random.randint(4)
                if np.random.uniform() > 0.5:
                    #signal
                    true = np.random.normal(1.0, 0.25)
                    syst = signalsyst
                else:
                    #bkgd
                    true = np.random.uniform(0.0, 5.0)
                    syst = bkgdsyst
                reco = np.random.normal(true, 0.0001)
                if reco > 0.0 and true > 0.0:
                    yield (reco, true, nupdg), 1.0, syst
        def gensk(N):
            eff = 0.5
            for (reco, true, nupdg), weight, syst in gen(N):
                   yield (reco, true, nupdg), eff*weight, weight, syst
        iternd280 = gen(10**5)
        itersuperk = gensk(1000)
        systematics = [("signal", [-5.0, 0.0, 5.0]), ("bkgd", [-5.0, 0.0, 5.0])]
        observables = ["reco_energy"]
        nd280 = BinnedSample("nd280", binning, observables, iternd280, cache_name="nd280_" + cachestr, systematics=systematics)
        prior = GaussianPrior(["signal", "bkgd"], [0.0, 0.0], [0.1, 0.1])
        if withosc:
            superk = BinnedSampleWithOscillation("superk", binning, observables, itersuperk, "true_energy", "true_nupdg", 295.0, cache_name="superk_" + cachestr)
            model = CombinedBinnedSample([nd280, superk])
            prior = CombinedPrior([OscillationParametersPrior(), prior])
        else:
            model = CombinedBinnedSample([nd280])
        #generate monte carlo
        ratevector = model
        toymc = ToyMC(ratevector, prior.generator)
        lhd_data = EventRateLikelihood(model, data=toymc.asimov().vec)
        lhd_prior = prior.likelihood
        lhd = SumLikelihood([lhd_data, lhd_prior])
        return model, toymc, lhd

def main():
    return unittest.main()

if __name__ == "__main__":
    main()

