"""Tools and definitions of colour palettes for ROOT.
"""

import ROOT
import numpy
import itertools

class ColorPalette:
    """Provides iterators over a set of colours for plot formatting.
    
    Defines various color schemes and provides functions to return iterable over 
    them. For example, if you have a set of 4 histograms, you may choose the 
    colour scheme "mode_redBlue". For example,
    
    >>> from commonAnalysis.plot import palette
    >>> 
    >>> colors = palette.ColorPalette.createPaletteCycle(palette.ColorPalette.mode_redBlue)
    >>> for hist,col in zip(listOfHistograms,colors):
    >>>     hist.SetLineColor(col)
    >>>
    """
    
    #enumerations for palettes
    mode_rootFirst = 1
    mode_redBlue = 2
    mode_violetGreenOrange = 3
    mode_tetrad = 4
    mode_cbSafe = 5
    #define default behaviour
    mode_default = 4
    
    #red and blue RGB numbers
    rgb_redBlue = [ (202, 0, 32), 
                   (244, 165, 130),
                   #(247, 247, 247), 
                   (146, 197, 222),
                   (5, 113, 176) 
                   ] 
    rgb_violetGreenOrange = [(250,187,0),
                             (96,128,0), 
                             (96,0,128),
                              (189,142,0),
                              (128,32,0),
                              (0,96,128),
                             ]
    
    rgb_tetrad = [ (0xFD,0x2B,0x0A),
                  (0xFD,0xB6,0x0A),
                  (0x1B,0x35,0xAD),
                  (0x07,0xBB,0x3E),
                  (0xA4,0x19,0x03),
                  (0xA4,0x75,0x03),
                  (0x09,0x1B,0x70),
                  (0x02,0x7A,0x26),
                             ]
    
    rgb_cbSafe = [ 
                   (0xFF,0xDA,0x00),
                   (0x42,0x12,0xAF),
                   (0xFF,0x74,0x00),
                   (0x00,0x99,0x99),
                   (0xA6,0x8E,0x00),
                   (0x27,0x06,0x72),
                   (0xA6,0x4B,0x00),
                   (0x00,0x63,0x63),
                  ]
    
    @staticmethod
    def createPaletteCycle(mode=mode_default):
        """Get iterator over palette colors.
        
        Once you reach the end of the iterator it will start again from the beginner i.e. it will iterate infinately. 
        
        :param mode: one of the modes defined in palette.ColorPalette
        :type mode: int
        :returns: returns an iterator over the colors in the requested palette.
        
        """
        return iter(ColorPalette(mode))
    
    def __init__(self, mode=mode_default):
        self.colorObjects = dict()
        self.colorCycle = itertools.cycle([1]) 
        self._createColors(mode)
    def __iter__(self):
        return self.colorCycle
    
    def _createColors(self, mode):
        cols = None
        if mode == ColorPalette.mode_rootFirst:
            cols = range(2,10) + [11,15,28,38,46,40,49]
        elif mode == ColorPalette.mode_redBlue:
            cols = self._convertRgbList(ColorPalette.rgb_redBlue)
        elif mode == ColorPalette.mode_violetGreenOrange:
            cols = self._convertRgbList(ColorPalette.rgb_violetGreenOrange)
        elif mode == ColorPalette.mode_tetrad:
            cols = self._convertRgbList(ColorPalette.rgb_tetrad)
        elif mode == ColorPalette.mode_cbSafe:
            cols = self._convertRgbList(ColorPalette.rgb_cbSafe)
        self.colorList = cols
        self.colorCycle = itertools.cycle(cols)
        
    def _convertRgbList(self, listOfRgb):
        result = []
        for rgbVal in listOfRgb:
            num = self._makeColor(*rgbVal)
            result.append(num)
        return result
    
    def _makeColor(self, r, g, b):
        #convert from hex (0-255) to ROOT (0.0-1.0)
        r = float(r)/255.0
        g = float(g)/255.0
        b = float(b)/255.0
        ROOT.TColor.InitializeColors()
        #get number for TColor
        #Search for the highest color index not used in ROOT:
        #We do not want to overwrite some colors...
        colorTable = ROOT.gROOT.GetListOfColors()
        lastColor = colorTable.Last()
        highestIndex = lastColor.GetNumber()
        color = lastColor
        while color:
             color = colorTable.Before(color)
             if color and color.GetNumber() > highestIndex:
                highestIndex = color.GetNumber();
        #create the color
        rgb = (.123, .456, .789)
        key = (r,g,b)
        col = ROOT.TColor(highestIndex+1,r,g,b,str(key))
        ROOT.SetOwnership(col,False)
        self.colorObjects[key] = col
        num = col.GetNumber()
        #print r,g,b,num,highestIndex
        return num
    
    def getColorList(self):
        return self.colorList

###############################################################################

class ColorGradient:
    """Used to set the ROOT color gradient used in 2D histograms drawn with COLZ.
    
    For example if you want a version that goes from Red-to-White-to-Blue 
    gradient (useful for the covariance matrices that T2K loves so much),
    
    >>> from commonAnalysis.plot import palette
    >>> palette.ColorGradient.setRedWhiteBlue()
    >>> #now draw and save plots
    >>>
      
    """
    @staticmethod
    def setPalette(red, green, blue, num=19):
        r = numpy.array(red)
        g = numpy.array(green)
        b = numpy.array(blue)
        nStops = len(red)
        uniform = [ float(i)/float(nStops-1) for i in xrange(nStops) ]
        stop = numpy.array(uniform)
        result = numpy.array(ROOT.TColor.CreateGradientColorTable(len(stop), stop, r, g, b, num))
        return
    
    @staticmethod
    def setDefault():
        ColorGradient.setInverseDarkBodyRadiator()
    
    @staticmethod
    def setRainbow():
        """Set the standard ROOT rainbow palette.""" 
        ROOT.gStyle.SetPalette(1)
    
    @staticmethod
    def setDarkBodyRadiator():
        """Set the standard ROOT Dark Body Radiator palette.""" 
        red = [ 0.00, 0.50, 1.00, 1.00, 1.00 ]
        green = [ 0.00, 0.00, 0.55, 1.00, 1.00 ]
        blue  = [ 0.00, 0.00, 0.00, 0.00, 1.00 ]
        ColorGradient.setPalette(red,green,blue)
        return
    
    @staticmethod
    def setInverseDarkBodyRadiator():
        """Set the standard ROOT Dark Body Radiator palette.""" 
        red = [ 0.00, 0.50, 1.00, 1.00, 1.00 ][::-1]
        green = [ 0.00, 0.00, 0.55, 1.00, 1.00 ][::-1]
        blue  = [ 0.00, 0.00, 0.00, 0.00, 1.00 ][::-1]
        ColorGradient.setPalette(red,green,blue)
        return
    
    @staticmethod
    def setRedWhiteBlue():
        """Palette that goes from (top) Red-White-Blue (bottom).
        Useful for covariance matrices.
        """
        r = [0.3, 1.0, 1.0]
        g = [0.3, 1.0, 0.3]
        b = [1., 1.0, 0.3]
        ColorGradient.setPalette(r,g,b)
        return
    
    @staticmethod
    def setRedBlue():
        """Palette that goes from (top) Red-Blue (bottom).
        """
        r = [0.0, 1.0]
        g = [0.0, 0.0]
        b = [1.0, 0.0]
        ColorGradient.setPalette(r,g,b)
        return
    
    @staticmethod
    def setGrayScale():
        """Palette that goes from (top)Black-White(bottom)"""
        r = [1.0, 0.4]
        g = [1.0, 0.4]
        b = [1.0, 0.4]
        ColorGradient.setPalette(r,g,b)
        return

###############################################################################

def _generate2DHistogram():
    hist = ROOT.TH2F("hist","hist",5,0,5,5,0,5)
    for ir in xrange(5):
        for ic in xrange(5):
            v = float(ir+1)+float(ic+1)/10.0
            hist.SetBinContent(ir+1,ic+1,v-2.5)
    return hist

def _test_ColorGradient():
    canv = ROOT.TCanvas()
    hist = _generate2DHistogram()
    ROOT.gStyle.SetOptStat(0)
    ColorGradient.setRedWhiteBlue()
    hist.Draw("COLZ,TEXT")
    canv.SetGrayscale()
    canv.Update()
    raw_input("wait")
    return

def _test_ColorPalette():
    import numpy
    canv = ROOT.TCanvas()
    cp = iter(ColorPalette(ColorPalette.mode_default))
    cols = [ cp.next() for i in range(8) ]
    print cols
    boxes = []
    hist = ROOT.TH1F("test","test",100,-3,3)
    hist.SetLineWidth(2)
    hist.Draw()
    for i,c in enumerate(cols):
        yMin = float(i)/float(len(cols))
        yMax = float(i+1)/float(len(cols))
        hist.FillRandom("gaus",10000)
        hist.SetLineColor(c)
        cHist = hist.Clone()
        cHist.Draw("SAME")
        boxes.append(cHist)
    #canv.Modified()
    #canv.Update()
    canv.SaveAs("blah.eps")
    canv.SaveAs("blah.png")
    raw_input("wait")
    return

###############################################################################

def main():
    _test_ColorGradient()
    #_test_ColorPalette()

if __name__ == "__main__":
    main()