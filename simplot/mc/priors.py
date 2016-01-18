import math

from simplot.mc.generators import GaussianGenerator, MultiVariateGaussianGenerator, UniformGenerator, GeneratorList
from simplot.mc.likelihood import GaussianLikelihood, MultiVariateGaussianLikelihood, CombinedLikelihood, ConstantLikelihood
from simplot.pdg import PdgNeutrinoOscillationParameters
from simplot.binnedmodel.sample import OscParMode

class Prior(object):
    def __init__(self, gen, lhd):
        self.generator = gen
        self.likelihood = lhd

class GaussianPrior(Prior):
    def __init__(self, parameter_names, mu, sigma, seed=None):
        gen = GaussianGenerator(parameter_names, mu, sigma, seed=seed)
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        super(GaussianPrior, self).__init__(gen, lhd)

class MultiVariateGaussianPrior(Prior):
    def __init__(self, parameter_names, mu, cov, seed=None):
        gen = MultiVariateGaussianGenerator(parameter_names, mu, cov, seed=seed)
        lhd = MultiVariateGaussianLikelihood(parameter_names, mu, cov)
        super(MultiVariateGaussianPrior, self).__init__(gen, lhd)

class UniformPrior(Prior):
    def __init__(self, parameter_names, value, range_, seed=None):
        gen = UniformGenerator(parameter_names, value, range_, seed=seed)
        lhd = ConstantLikelihood(parameter_names, 0.0)
        super(UniformPrior, self).__init__(gen, lhd)

class CombinedPrior(Prior):
    def __init__(self, priors):
        gen = [p.generator for p in priors]
        gen = GeneratorList(*gen)
        lhd = CombinedLikelihood([p.likelihood for p in priors])
        super(CombinedPrior, self).__init__(gen, lhd)

class OscillationParametersPrior(CombinedPrior):
    def __init__(self, values=None, usereactorconstraint=False, seed=None, oscparmode=OscParMode.SINSQTHETA):
        if values is None:
            values = PdgNeutrinoOscillationParameters()
        priors = self._buildpriors(values, usereactorconstraint=usereactorconstraint, seed=seed, oscparmode=oscparmode)
        super(OscillationParametersPrior, self).__init__(priors)
        
    def _buildpriors(self, values, usereactorconstraint, seed=None, oscparmode=OscParMode.SINSQTHETA):
        if seed is not None:
            seed = seed + 230932904 # offset seed, incase user is giving sequential seeds to various generators.
        priors = []
        if oscparmode == OscParMode.SINSQTHETA:
            parnames = type(values).ALL_PARS_SINSQ
        elif oscparmode == OscParMode.SINSQ2THETA:
            parnames = type(values).ALL_PARS_SINSQ2
        else:
            parnames = type(values).ALL_PARS
        for p in parnames:
            if seed is not None:
                seed += 1
            val = values.value(p)
            err = values.error(p)
            if p == type(values).DELTACP:
                p = UniformPrior([p], [val], range_=[(-math.pi, math.pi)], seed=seed)
            elif p == type(values).THETA13 and not usereactorconstraint:
                p = UniformPrior([p], [val], range_=[(-math.pi, math.pi)], seed=seed)
            elif (p == type(values).SINSQ2THETA13 or p == type(values).SINSQTHETA13) and not usereactorconstraint:
                p = UniformPrior([p], [val], range_=[(0.0, 1.0)], seed=seed)
            else:
                p = GaussianPrior([p], mu=[val], sigma=[err], seed=seed)
            priors.append(p)
        return priors
