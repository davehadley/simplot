""" plot package provides tools for generating, drawing and saving histograms.

"""

import ROOT#@UnresolvedImport
ROOT.PyConfig.IgnoreCommandLineOptions = True # Prevent ROOT from hi-jacking --help

import drawoptions
import drawtools
import rootio
import ntuple
import palette
import style
