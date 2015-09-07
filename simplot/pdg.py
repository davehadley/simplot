import math


class PdgNeutrinoOscillationParameters(object):

    THETA12 = "theta12"
    THETA23 = "theta23"
    THETA13 = "theta13"
    DELTACP = "deltacp"
    SDM = "sdm"
    LDM = "ldm"

    SINSQ2THETA12 = "sinsq2theta12"
    SINSQ2THETA23 = "sinsq2theta23"
    SINSQ2THETA13 = "sinsq2theta13"

    ALL_PARS = [THETA12, THETA23, THETA13, DELTACP, SDM, LDM]
    ALL_PARS_SINSQ2 = [SINSQ2THETA12, SINSQ2THETA23, SINSQ2THETA13, DELTACP, SDM, LDM]

    NUM_OSC_PARS = 6

    def __init__(self):
        num_osc_pars = 6
        #osc values
        sinsq_2theta23 = 1.0
        sinsq_2theta13 = 0.095
        sinsq_2theta12 = 0.857
        deltacp = 0.0 #math.pi / 2.0
        delta_msq32 = 2.32e-3
        delta_msq12 = 7.5e-5
        
        #osc errors
        # PDG error on theta23 gives sinsq_2theta23 > 0.95 at 90% confidence.
        # When modelling with a half-Gaussian this limit corresponds to:
        #  Gaus(x; mu=1.0, sigma=0.03 | x < 1.0)
        sinsq_2theta23_errors = 0.03
        sinsq_2theta13_errors = 0.010
        sinsq_2theta12_errors = 0.024
        deltacp_errors = 0.0
        delta_msq32_errors = 0.10e-3
        delta_msq12_errors = 0.2e-5
        
        #transformed
        theta12 = math.asin(math.sqrt(sinsq_2theta12)) / 2.0 #in radians
        theta13 = math.asin(math.sqrt(sinsq_2theta13)) / 2.0 #in radians
        theta23 = math.asin(math.sqrt(sinsq_2theta23)) / 2.0 #in radians
        sdm = delta_msq12 #in eV^2
        ldm = delta_msq32 #in eV^2
        #store values
        cls = type(self)
        self._values = {cls.THETA12: theta12,
                        cls.THETA23: theta13,
                        cls.THETA13: theta23,
                        cls.SINSQ2THETA12: sinsq_2theta12,
                        cls.SINSQ2THETA23: sinsq_2theta23,
                        cls.SINSQ2THETA13: sinsq_2theta13,
                        cls.DELTACP: 0.0,
                        cls.LDM: ldm,
                        cls.SDM: sdm,
        }
        self._errors = {cls.SINSQ2THETA12: sinsq_2theta12_errors,
                        cls.SINSQ2THETA23: sinsq_2theta23_errors,
                        cls.SINSQ2THETA13: sinsq_2theta13_errors,
                        cls.DELTACP: 0.0,
                        cls.LDM: ldm,
                        cls.SDM: sdm,
        }

    def value(self, par):
        return self._values[par]

    def error(self, par):
        return self._errors[par]        



