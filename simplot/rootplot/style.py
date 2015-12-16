# -*- coding: utf-8 -*-
"""Defines several ROOT styles.

Two main styles are provided: ATLAS style, T2KStyle. 
ATLAS style is recommended (it is much better than T2K style).
Use this once at the start of your plotting script. For example,

>>> from commonAnalysis.plot import style
>>> style.setATLASStyle()

"""

import ROOT

###############################################################################

def setT2KStyle():
    print "setting T2K style"
    canv = ROOT.TCanvas("temp","temp",10,10)
    style = ROOT.gROOT.GetStyle("T2KStyle")
    if not style:
        #create new style
        style = ROOT.TStyle("T2KStyle","T2K style.")
        ROOT.SetOwnership(style,False)
        style = configOfficialATLASStyle(style)
        #my modifications
        style.SetPalette(1) #rainbow palette 
        #style.SetLegendBorderSize(0);#make legend without border
        style.SetPadGridX(0);
        style.SetPadGridY(0);

        #Canvas Options
        style.SetCanvasBorderSize(0);
        #width = 1200;
        #style.SetCanvasDefW(width)
        #style.SetCanvasDefH(int(double(width)/1.618));#golden ratio ~1.618:1
        #style.SetCanvasDefH(int(float(width)*(3./4.))) #4:3 ratio

        #Title Options
        style.SetTitleFont(42)
        style.SetTitleBorderSize(0);#no title border
        style.SetTitleFillColor(0);
        style.SetStatFont(42)
        style.SetStatBorderSize(0)
        style.SetPadTopMargin(0.10);
        style.SetOptStat(0) #style.SetOptStat("eou"); #style.SetOptStat(10);#. iourmen n=name e=entries m=mean r=rms u=underflow o=overflow i=integral
        
        #testing
        style = configOfficialT2KStyle(style)
        
#        #overwrite things that need overwriting
        style.SetPadLeftMargin(0.16);
#        tsize=0.05;
#        style.SetTextSize(tsize);
#        style.SetLabelSize(tsize,"x");
#        style.SetTitleSize(tsize,"x");
#        style.SetLabelSize(tsize,"y");
#        style.SetTitleSize(tsize,"y");
#        style.SetLabelSize(tsize,"z");
#        style.SetTitleSize(tsize,"z");

        textScale = 1.1
        #. use large Times-Roman fonts
        #style.SetTextFont(132);
        style.SetTextSize(0.08*textScale);
        style.SetLabelFont(132,"x");
        style.SetLabelFont(132,"y");
        style.SetLabelFont(132,"z");
        style.SetLabelSize(0.05*textScale,"x");
        style.SetTitleSize(0.06*textScale,"x");
        style.SetLabelSize(0.05*textScale,"y");
        style.SetTitleSize(0.06*textScale,"y");
        style.SetLabelSize(0.05*textScale,"z");
        style.SetTitleSize(0.06*textScale,"z");
        #style.SetLabelFont(132,"t");
        #style.SetTitleFont(132,"x");
        #style.SetTitleFont(132,"y");
        #style.SetTitleFont(132,"z");
        #style.SetTitleFont(132,"t"); 
        style.SetTitleFillColor(0);
        style.SetTitleX(0.25*textScale);
        style.SetTitleFontSize(0.08*textScale);
        #style.SetTitleFont(132,"pad");
                
        style.SetLegendFont(132);
        
        style.SetTitleXOffset(0.8*style.GetTitleXOffset());
        style.SetTitleYOffset(0.8*style.GetTitleYOffset());
        
        #lines
        style.SetLineScalePS(4.);
        style.SetHistLineWidth(4);
        
    style.cd()
    ROOT.gROOT.ForceStyle() 

###############################################################################

def setOfficialT2KStyle():
    print "setting T2K style"
    canv = ROOT.TCanvas("temp","temp",10,10)
    style = ROOT.gROOT.GetStyle("T2KOfficialStyle")
    if not style:
        style = ROOT.TStyle("T2KOfficialStyle","Official T2K style.")
        ROOT.SetOwnership(style,False)
        style = configOfficialT2KStyle(style)
    style.cd()
    ROOT.gROOT.ForceStyle()

###############################################################################

def setATLASStyle():
    #print "setting ATLAS style"
    canv = ROOT.TCanvas("temp","temp",10,10)
    style = ROOT.gROOT.GetStyle("AtlasStyle")
    if not style:
        #create new style
        style = ROOT.TStyle("AtlasStyle","ZbbNtuple style based on ATLAS style")
        style = configOfficialATLASStyle(style)
        #my modifications
        style.SetPalette(1) #rainbow palette 
        #style.SetFillStyle(0);#.make added objects transparent
        style.SetLegendBorderSize(1);#make legend without border
        style.SetPadGridX(0);
        style.SetPadGridY(0);

        #Canvas Options
        style.SetCanvasBorderSize(0);
        #width = 1900;
        #style.SetCanvasDefW(width);#.fixed width 1200
        #style.SetCanvasDefH(int(double(width)/1.618));#.golden ratio ~1.618:1

        #Title Options
        #style.SetOptTitle(0);#show title
        style.SetTitleBorderSize(0);#no title border
        style.SetTitleFillColor(0);
        style.SetTitleStyle(0);
        #style.SetTitleX( 1.0 );
        #style.SetTitleY( 1.0 );
        style.SetTitleFont(42);
        #style.SetPadTopMargin(0.20);
        #Stats box Options
        #style.SetStatFont(42);
        #style.SetStatFontSize(0.04);
        #style.SetOptStat(0); #style.SetOptStat(10);#. iourmen n=name e=entries m=mean r=rms u=underflow o=overflow i=integral
        style.SetOptStat(0)#style.SetOptStat("eou")
        
        #Markers
        #style.SetErrorX(0.0);#.Removes horizontal bin errors
        #style.SetMarkerSize(0.5);

        #lines
        style.SetLineScalePS(4.);
        style.SetHistLineWidth(4);

        style.cd();
        ROOT.SetOwnership(style,False)
    style.cd()
    ROOT.gROOT.ForceStyle()
    
###############################################################################

def configOfficialATLASStyle(atlasStyle):
    icol = 0 #color = white
    atlasStyle.SetFrameBorderMode(icol);
    atlasStyle.SetFrameFillColor(icol);
    atlasStyle.SetCanvasBorderMode(icol);
    atlasStyle.SetCanvasColor(icol);
    atlasStyle.SetPadBorderMode(icol);
    atlasStyle.SetPadColor(icol);
    atlasStyle.SetStatColor(icol);
    atlasStyle.SetFillColor(icol); # don't use: white fill color for *all* objects
    
    # set the paper size
    atlasStyle.SetPaperSize(20,26);
    # set margin sizes
    atlasStyle.SetPadTopMargin(0.05);
    atlasStyle.SetPadRightMargin(0.05);
    atlasStyle.SetPadBottomMargin(0.16);
    atlasStyle.SetPadLeftMargin(0.16);
    
    # set title offsets (for axis label)
    atlasStyle.SetTitleXOffset(1.4);
    atlasStyle.SetTitleYOffset(1.4);
    
    # use large fonts
    #Int_t font=72; # Helvetica italics
    font=42; # Helvetica
    tsize=0.05;
    atlasStyle.SetTextFont(font);
    atlasStyle.SetTextSize(tsize);
    atlasStyle.SetLabelFont(font,"x");
    atlasStyle.SetTitleFont(font,"x");
    atlasStyle.SetLabelFont(font,"y");
    atlasStyle.SetTitleFont(font,"y");
    atlasStyle.SetLabelFont(font,"z");
    atlasStyle.SetTitleFont(font,"z");
    atlasStyle.SetLabelSize(tsize,"x");
    atlasStyle.SetTitleSize(tsize,"x");
    atlasStyle.SetLabelSize(tsize,"y");
    atlasStyle.SetTitleSize(tsize,"y");
    atlasStyle.SetLabelSize(tsize,"z");
    atlasStyle.SetTitleSize(tsize,"z");
    
    # use bold lines and markers
    atlasStyle.SetMarkerStyle(20);
    atlasStyle.SetMarkerSize(1.2);
    atlasStyle.SetHistLineWidth(2);
    atlasStyle.SetLineStyleString(2,"[12 12]"); # postscript dashes

    # get rid of X error bars
    #atlasStyle.SetErrorX(0.001);
    # get rid of error bar caps
    atlasStyle.SetEndErrorSize(0.);

    # do not display any of the standard histogram decorations
    atlasStyle.SetOptTitle(0);
    #atlasStyle.SetOptStat(1111);
    atlasStyle.SetOptStat(0);
    #atlasStyle.SetOptFit(1111);
    atlasStyle.SetOptFit(0);

    #put tick marks on top and RHS of plots
    atlasStyle.SetPadTickX(1);
    atlasStyle.SetPadTickY(1);

    return atlasStyle;

###############################################################################

def configOfficialT2KStyle(t2kStyle):
    #T2K style definition
    #Adopted from BaBar collaboration
    #Add the following lines to the start of your rootlogon.C file
    #TStyle *t2kStyle= new TStyle("T2K","T2K approved plots style");

    #use plain black on white colors
    t2kStyle.SetFrameBorderMode(0);
    t2kStyle.SetCanvasBorderMode(0);
    t2kStyle.SetPadBorderMode(0);
    t2kStyle.SetPadColor(0);
    t2kStyle.SetCanvasColor(0);
    t2kStyle.SetStatColor(0);
    t2kStyle.SetFillColor(0);
    t2kStyle.SetLegendBorderSize(1); 
    
    #. set the paper & margin sizes
    t2kStyle.SetPaperSize(20,26);
    t2kStyle.SetPadTopMargin(0.05);
    t2kStyle.SetPadRightMargin(0.05);
    t2kStyle.SetPadBottomMargin(0.16);
    t2kStyle.SetPadLeftMargin(0.12);
    
    #. use large Times-Roman fonts
    t2kStyle.SetTextFont(132);
    t2kStyle.SetTextSize(0.08);
    t2kStyle.SetLabelFont(132,"x");
    t2kStyle.SetLabelFont(132,"y");
    t2kStyle.SetLabelFont(132,"z");
    t2kStyle.SetLabelSize(0.05,"x");
    t2kStyle.SetTitleSize(0.06,"x");
    t2kStyle.SetLabelSize(0.05,"y");
    t2kStyle.SetTitleSize(0.06,"y");
    t2kStyle.SetLabelSize(0.05,"z");
    t2kStyle.SetTitleSize(0.06,"z");
    t2kStyle.SetLabelFont(132,"t");
    t2kStyle.SetTitleFont(132,"x");
    t2kStyle.SetTitleFont(132,"y");
    t2kStyle.SetTitleFont(132,"z");
    t2kStyle.SetTitleFont(132,"t"); 
    t2kStyle.SetTitleFillColor(0);
    t2kStyle.SetTitleX(0.25);
    t2kStyle.SetTitleFontSize(0.08);
    t2kStyle.SetTitleFont(132,"pad");
    
    #. use bold lines and markers
    t2kStyle.SetMarkerStyle(20);
    #t2kStyle.SetHistLineWidth(1.85);
    t2kStyle.SetHistLineWidth(2);#increased to 2 to make it work with python.
    t2kStyle.SetLineStyleString(2,"[12 12]"); #. postscript dashes
    
    #. get rid of X error bars and y error bar caps
    t2kStyle.SetErrorX(0.001);
    
    
    
    #. do not display any of the standard histogram decorations
    t2kStyle.SetOptTitle(0);
    t2kStyle.SetOptStat(0);
    t2kStyle.SetOptFit(0);
    
    #. put tick marks on top and RHS of plots
    t2kStyle.SetPadTickX(1);
    t2kStyle.SetPadTickY(1);
    
    #TODO: make python compatible
    #. Add a greyscale palette for 2D plots
    #ncol=50;
    #dcol = 1./float(ncol);
    #gray = 1;
    #TColor **theCols = new TColor*[ncol];
    #for (int i=0;i<ncol;i++) theCols[i] = new TColor(999-i,0.0,0.7,0.7);
    #for (int j = 0; j < ncol; j++) {
    #  theCols[j].SetRGB(gray,gray,gray);
    #  gray -= dcol;
    #}
    #int ColJul[100];
    #for  (int i=0; i<100; i++) ColJul[i]=999-i;
    #t2kStyle.SetPalette(ncol,ColJul);
    
    
    
    #. Define a nicer color palette (red.blue)
    #. Uncomment these lines for a color palette (default is B&W)
    #. t2kStyle.SetPalette(1,0);  #. use the nice red.blue palette
    #. const Int_t NRGBs = 5;
    #. const Int_t NCont = 255;
    #.
    #. Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    #. Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    #. Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    #. Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    #. TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
    #. NCont);
    #. t2kStyle.SetNumberContours(NCont); 
    
    #. End of definition of t2kStyle
    return t2kStyle

###############################################################################

def setECalPidStyle():
    setATLASStyle()
    
###############################################################################

def setDefaultStyle():
    setATLASStyle()
    
###############################################################################

#set default style
#setDefaultStyle()

    
