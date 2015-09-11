from simplot.binnedmodel.xsecweights import XsecWeights, InterpolatedWeightCalc
from simplot.binnedmodel.fluxweights import FluxWeights

################################################################################

class Systematics(object):
    @property
    def spline_parameter_values(self):
        return []

    @property
    def parameter_names(self):
        raise NotImplementedError("ERROR: child class must implement this method.")

    def __call__(self, parameter_names, systhist, nominalhist):
        raise NotImplementedError("ERROR: child class must implement this method.")        

################################################################################

class SplineSystematics(Systematics):

    def __init__(self, spline_parameter_values):
        self._spline_parameter_values = spline_parameter_values

    def __call__(self, parameter_names, systhist, nominalhist):
        xsec_weights = self._buildxsecweights(self.spline_parameter_values, systhist, nominalhist)
        flux_weights = self._buildfluxweights()
        return xsec_weights, flux_weights

    @property
    def spline_parameter_values(self):
        return self._spline_parameter_values

    @property
    def parameter_names(self):
        return self._build_parameter_names(self._spline_parameter_values, None)

    def _buildxsecweights(self, systematicsvalues, systhist, hist):
        xsecweights = None
        if systematicsvalues is not None:
            wclist = []
            for isyst, (syst, parval) in enumerate(systematicsvalues):
                assert len(parval) == len(systhist[isyst])
                l = zip(parval, systhist[isyst])
                l.sort() # sort by parameter value
                weights = [x[1].array() for x in l] 
                parval = [x[0] for x in l]
                wc = InterpolatedWeightCalc(hist.array(), parval, weights, syst, self.parameter_names)
                wclist.append(wc)
            xsecweights = XsecWeights(hist.array(), wclist)
        return xsecweights

    def _buildfluxweights(self):
        fw = None
        return fw

    def _build_parameter_names(self, systematics, fluxsystematics):
        parameter_names = []
        if systematics:
            for s, _ in systematics:
                parameter_names.append(s)
        if fluxsystematics:
            for s in fluxsystematics.parameter_names:
                parameter_names.append(s)
        return parameter_names

################################################################################

class FluxSystematics(Systematics):
    def __init__(self, enudim, nupdgdim, fluxparametermap):
        self._dim_enutrue = enudim
        self._dim_nupdg = nupdgdim
        self._fluxparametermap = fluxparametermap

    @property
    def parameter_names(self):
        return self._fluxparametermap.keys()

    def __call__(self, parameter_names, systhist, nominalhist):
        xsec_weights = None
        flux_weights = self._buildfluxweights(nominalhist)
        return xsec_weights, flux_weights

    def _buildfluxweights(self, nominalhist):
        fw = FluxWeights(parameter_names, nominalhist.array().shape(),
                     self._dim_enutrue, 
                     self._dim_nupdg,
                     -1, # no detector dimension
                     binning.dim_beammode,
                     parametermap=self._fluxparametermap,
        )
        return

################################################################################

class FluxAndSplineSystematics(Systematics):
    def __init__(self, spline_parameter_values):
        self._splinesyst = SplineSystematics(spline_parameter_values)
        self._fluxsyst = FluxSystematics()

    def __call__(self, parameter_names, systhist, nominalhist):
        xsec_weights, _ = self._splinesyst(parameter_names, systhist, nominalhist)
        _, flux_weights = self._fluxsyst(parameter_names, systhist, nominalhist)

    @property
    def parameter_names(self):
        return self._splinesyst.parameter_names + self._fluxsyst.parameter_names

################################################################################
