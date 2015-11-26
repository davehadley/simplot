from collections import OrderedDict

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
        xsec_weights = self._buildxsecweights(self.spline_parameter_values, parameter_names, systhist, nominalhist)
        flux_weights = self._buildfluxweights()
        det_weights = self._builddetweights()
        return det_weights, xsec_weights, flux_weights

    @property
    def spline_parameter_values(self):
        return self._spline_parameter_values

    @property
    def parameter_names(self):
        return self._build_parameter_names(self._spline_parameter_values)

    def _buildxsecweights(self, systematicsvalues, parameter_names, systhist, hist):
        xsecweights = None
        if systematicsvalues is not None:
            wclist = []
            for isyst, (syst, parval) in enumerate(systematicsvalues):
                assert len(parval) == len(systhist[isyst])
                l = zip(parval, systhist[isyst])
                l.sort() # sort by parameter value
                weights = [x[1].array() for x in l] 
                parval = [x[0] for x in l]
                wc = InterpolatedWeightCalc(hist.array(), parval, weights, syst, parameter_names)
                wclist.append(wc)
            xsecweights = XsecWeights(hist.array(), wclist)
        return xsecweights

    def _buildfluxweights(self):
        fw = None
        return fw

    def _builddetweights(self):
        dw = None
        return dw

    def _build_parameter_names(self, systematics):
        parameter_names = []
        if systematics:
            for s, _ in systematics:
                parameter_names.append(s)
        return parameter_names

################################################################################

class FluxSystematics(Systematics):
    def __init__(self, enudim, nupdgdim, beammodedim, fluxparametermap):
        self._dim_enutrue = enudim
        self._dim_nupdg = nupdgdim
        self._dim_beammode = beammodedim
        self._fluxparametermap = fluxparametermap

    @property
    def parameter_names(self):
        return self._fluxparametermap.keys()

    def __call__(self, parameter_names, systhist, nominalhist):
        det_weights = None
        xsec_weights = None
        flux_weights = self._buildfluxweights(parameter_names, nominalhist)
        return det_weights, xsec_weights, flux_weights

    def _buildfluxweights(self, parameter_names, nominalhist):
        fw = FluxWeights(parameter_names, nominalhist.array().shape(),
                     self._dim_enutrue, 
                     self._dim_nupdg,
                     -1, # no detector dimension
                     self._dim_beammode,
                     parametermap=self._fluxparametermap,
        )
        return fw

    @classmethod
    def make_flux_parameter_map(cls, enubinning, flux_error_binning, name_pattern=None):
        cls._check_bin_mapping(enubinning, flux_error_binning)
        result = OrderedDict()
        detbin = -1
        for key, beambin, flavbin, binning in flux_error_binning:
            for ipar, (low, high) in enumerate(zip(binning[:-1], binning[1:])):
                index = tuple(list(key) + [ipar])
                if name_pattern is None:
                    name_pattern = "f_" + "_".join(["%s"] * len(index))
                parname = name_pattern % index
                l = []
                for xi, xlo in enumerate(enubinning[:-1]):
                    if low <= xlo < high:
                        l.append((detbin, beambin, flavbin, xi))
                result[parname] = l
        return result

    @classmethod
    def _check_bin_mapping(cls, enubinning, flux_error_binning):
        for _, _, _, errbinning in flux_error_binning:
            for x in errbinning:
                if not x in enubinning:
                    raise Exception("ERROR: flux parameter mis-mapping. An error bin edge is not in the global enu binning", x, enubinning)
        return

################################################################################

class FluxAndSplineSystematics(Systematics):
    def __init__(self, spline_parameter_values, enudim, nupdgdim, beammodedim, fluxparametermap):
        self._splinesyst = SplineSystematics(spline_parameter_values)
        self._fluxsyst = FluxSystematics(enudim, nupdgdim, beammodedim, fluxparametermap)

    @property
    def spline_parameter_values(self):
        return self._splinesyst.spline_parameter_values

    def __call__(self, parameter_names, systhist, nominalhist):
        det_weights = None
        _, xsec_weights, _ = self._splinesyst(parameter_names, systhist, nominalhist)
        _, _, flux_weights = self._fluxsyst(parameter_names, systhist, nominalhist)
        return det_weights, xsec_weights, flux_weights

    @property
    def parameter_names(self):
        return self._splinesyst.parameter_names + self._fluxsyst.parameter_names

################################################################################
