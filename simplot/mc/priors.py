import math

from simplot.mc.generators import GaussianGenerator, UniformGenerator, GeneratorList
from simplot.mc.likelihood import GaussianLikelihood, CombinedLikelihood, ConstantLikelihood
from simplot.pdg import PdgNeutrinoOscillationParameters

class Prior(object):
    def __init__(self, gen, lhd):
        self.generator = gen
        self.likelihood = lhd

class GaussianPrior(Prior):
    def __init__(self, parameter_names, mu, sigma, seed=None):
        gen = GaussianGenerator(parameter_names, mu, sigma, seed=seed)
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        super(GaussianPrior, self).__init__(gen, lhd)

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
    def __init__(self, values=None, usereactorconstraint=False):
        if values is None:
            values = PdgNeutrinoOscillationParameters()
        priors = self._buildpriors(values, usereactorconstraint=usereactorconstraint)
        super(OscillationParametersPrior, self).__init__(priors)
        
    def _buildpriors(self, values, usereactorconstraint):
        priors = []
        parnames = type(values).ALL_PARS_SINSQ2
        for p in parnames:
            val = values.value(p)
            err = values.error(p)
            if p == type(values).DELTACP:
                p = UniformPrior([p], [val], range_=[(-math.pi, math.pi)])
            elif p == type(values).THETA13 and not usereactorconstraint:
                p = UniformPrior([p], [val], range_=[(-math.pi, math.pi)])
            elif p == type(values).SINSQ2THETA13 and not usereactorconstraint:
                p = UniformPrior([p], [val], range_=[(0.0, 1.0)])
            else:
                p = GaussianPrior([p], mu=[val], sigma=[err])
            priors.append(p)
        return priors
