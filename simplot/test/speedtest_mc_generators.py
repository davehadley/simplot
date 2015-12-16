from argparse import ArgumentParser
from simplot.profile import profile_func

import numpy as np

from simplot.mc.generators import GaussianGenerator, MultiVariateGaussianGenerator, ConstantGenerator, GeneratorList

_NPE = 10**6
_NPARS = 10

def speedtest_gaus(npe=_NPE, npars=_NPARS):
    n = [str(x) for x in range(npars)]
    mu = range(npars)
    sigma = range(npars)
    gen = GaussianGenerator(n, mu, sigma)
    for _ in xrange(npe):
       gen() 
    return

def speedtest_multigaus(npe=_NPE, npars=_NPARS):
    n = [str(x) for x in range(npars)]
    mu = range(npars)
    sigma = range(npars)
    cov = np.diag(sigma)
    gen = MultiVariateGaussianGenerator(n, mu, cov)
    for _ in xrange(npe):
       gen() 
    return

def speedtest_const(npe=_NPE, npars=_NPARS):
    n = [str(x) for x in range(npars)]
    mu = range(npars)
    gen = ConstantGenerator(n, mu)
    for _ in xrange(npe):
       gen() 
    return

def speedtest_list(npe=_NPE, npars=_NPARS):
    names = [str(x) for x in range(npars)]
    mu = range(npars)
    sigma = range(npars)
    gaus = GaussianGenerator(["gaus"+n for n in names], mu, sigma)
    cov = np.diag(sigma)
    multigaus = MultiVariateGaussianGenerator(["multi" + n for n in names], mu, cov)
    gen = GeneratorList(gaus, multigaus)
    for _ in xrange(npe):
       gen() 
    return


_GENERATOR_CHOICES = { "gaus" : speedtest_gaus,
                       "multigaus" : speedtest_multigaus,
                       "const" : speedtest_const,
                       "list" : speedtest_list,
}

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("generator", choices=_GENERATOR_CHOICES.keys(), default="gaus")
    return parser.parse_args()

def main():
    args = parsecml()
    func = _GENERATOR_CHOICES[args.generator]
    profile_func(func)
    return

if __name__ == "__main__":
    main()

