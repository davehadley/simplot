import ROOT
import StringIO
import pprint
import itertools
import os

import prettytable

###############################################################################

class StringFormatter:
    def __init__(self, fstr):
        self._fstr = fstr
    def __call__(self, entry):
        return self._fstr.format(entry)

_defaultformatter = StringFormatter("{:.2g}")

###############################################################################

class Table:
    def __init__(self, name, headrow):
        self._name = name
        self._headrow = None
        self._setheadrow(headrow)
        self._rows = []
            
    def _setheadrow(self, headrow):
        #ensure that header row is all strings
        if headrow is not None:
            headrow = [str(e) for e in headrow]
        self._headrow = headrow
        return
        
    def addrow(self, row):
        self._rows.append(row)
    
    def get_nrows(self):
        return len(self._rows)
    
    def get_ncols(self):
        return max( (len(l) for l in self._rows) )
    
    def _hasrowlabels(self):
        return any( (len(r)>0 and type(r[0]) is str for r in self._rows) )
        
    @staticmethod
    def tablefromdict(name, container, headrow=None):
        table = Table(name, headrow=headrow)
        if headrow:
            table._setheadrow(headrow)
        for key, row in container.iteritems():
            inputrow = [key] + row
            table.addrow(inputrow)
        return table
    
    @staticmethod
    def tablefrommatrix(name, matrix, xlabels=None, ylabels=None):
        nrows = matrix.GetNrows()
        ncols = matrix.GetNcols()
        if xlabels is None:
            xlabels = [str(i) for i in range(ncols)]
        if ylabels is None:
            ylabels = [str(i) for i in range(nrows)]
        tab = Table(name, headrow=ylabels)
        for j in xrange(nrows):
            row = [ xlabels[j] ]
            for i in xrange(ncols):
                v = matrix[i][j]
                row.append(v)
            tab.addrow(row)
        return tab

    def __str__(self):
        return "Table("+self.name+")"
    
###############################################################################

class TableOutputBase(object):
    def __init__(self, formatter):
        if formatter is None:
            formatter = _defaultformatter
        self._formatter = formatter
    
    def write(self, table, filename):
        try:
            os.makedirs(os.path.dirname(filename))
        except os.error:
            pass
        outfile = open(filename,"w")
        print >>outfile,self.getstring(table)
    
    def getstring(self, table):
        raise NotImplementedError("Users of this class should override this method.")

    def _get_printed_head_row(self, table):
        result = None
        if table._headrow:
            if len(table._headrow)==table.get_ncols()-1:
                result = [ "" ] + table._headrow
            else:
                result = list(table._headrow)
        return result

    def _convert_row_to_string(self, row):
        strRow = [self._convert_entry_to_string(v) for v in row]
        return strRow

    def _convert_entry_to_string(self, entry):
        strValue = None
        if isinstance(entry, basestring):
            #already a string, copy value
            strValue = str(entry)
        else:
            #not a string, try to get value and format
            try:
                value,error = self._get_value_and_error(entry)
                vs = self._formatter(value)
                es = ""
                if error is not None:
                    es = self._formatter(error)
                    strValue = vs + " +- " + es
                else:
                    strValue = vs
            except Exception:
                #can't convert to string, fall back on standard python string conversion
                strValue = str(entry)
                raise
        return strValue

    def _get_value_and_error(self, entry):
        value,error = None,None
        try:
            #is it a single number?
            value = float(entry)
        except:
            try:
                #try unpacking a tuple
                value, error = entry
                value = float(value)
                error = float(error)
            except:
                #don't know what to do now, raise an exception
                raise Exception("can't convert entry to value and error", entry)
        return value, error

###############################################################################

class TableLatexOutput(TableOutputBase):
    def __init__(self, formatter=None):
        super(TableLatexOutput, self).__init__(formatter)
        
    def getstring(self, table):
        return self._getlatex(table)
    
    def _getlatex(self, table):
        latex = StringIO.StringIO()
        nrows = table.get_nrows()
        ncols = table.get_ncols()
        newLine = " \\\\\n"
        colformat = "|".join(["c"]*ncols)
        print >>latex,"{"
        print >>latex,"%\\tiny" 
        print >>latex,"\\begin{tabular}{"+colformat+"}"
        hrow = self._get_printed_head_row(table)
        if hrow:
            allData = [hrow]+self.rows
        else:
            allData = table._rows
        for row in allData:
            strRow = self._convert_row_to_string(row)
            strRow = self._sanitise_latex_row(strRow)
            line = " & ".join(strRow) + newLine
            print >>latex,line,
        print >>latex,"\\end{tabular}"
        print >>latex,"}" 
        return latex.getvalue()
    
    def _sanitise_latex_row(self, row):
        result = [self._sanitise_latex_string(entry) for entry in row]
        return result
    
    def _sanitise_latex_string(self, entry):
        result = entry
        #deal with underscores
        result = result.replace("\\_","_")
        result = result.replace("_","\\_")
        #change \pm
        result = result.replace("+-","\\ensuremath{\\pm}")
        return result

###############################################################################

class TableAsciiOutput(TableOutputBase):
    def __init__(self, formatter=None):
        super(TableAsciiOutput, self).__init__(formatter)
        
    def getstring(self, table):
        return self._get_ascii_table(table)

    def _get_ascii_table(self, table):
        pt = prettytable.PrettyTable(self._get_printed_head_row(table))
        for row in table._rows:
            row = [ self._convert_entry_to_string(n) for n in row]
            pt.add_row(row)
        return str(pt)

###############################################################################

class TableHistogramOutput:
    def __init__(self, formatter=None):
        super(TableHistogramOutput, self).__init__(formatter)

    def drawColzPlot(self, minZ=None, maxZ=None, textFormat=None):
        name = self.name
        canv = ROOT.TCanvas(name,name,800,600)
        nRows = self.getNRows()
        nCols = self.getNCols()
        hasRowLabels = self.hasRowLabels()
        nY = nRows
        nX = nCols
        if hasRowLabels:
            nX = nCols - 1
        hist = ROOT.TH2F("hist"+name,"",nX,0,nX,nY,0,nY)
        hist.SetDirectory(0)
        autoMax = minZ
        for iX, iY in itertools.product(range(nX),range(nY)):
            if hasRowLabels:
                indexX = iX + 1
            else:
                indexX = iX
            indexY = iY
            binX = iX + 1
            binY = nRows - indexY
            entry = self.rows[indexY][indexX]
            #convert to floating point if we can
            value,error = self._getValueAndError(entry)
            autoMax = max(value,autoMax)
            hist.SetBinContent(binX,binY,value)
            if error is not None:
                hist.SetBinError(binX,binY,value)
        #set axis labels
        for i,xLabel in enumerate(self.headerRow):
            hist.GetXaxis().SetBinLabel(i+1,xLabel)
        if hasRowLabels:
            for iY in xrange(nY):
                indexY = iY
                yLabel = self.rows[indexY][0]
                binY = nRows - indexY
                hist.GetYaxis().SetBinLabel(binY,yLabel)
        if textFormat:
            startingTextFormat = ROOT.gStyle.GetPaintTextFormat()
            ROOT.gStyle.SetPaintTextFormat(textFormat)
        if minZ is not None:
            hist.SetMaximum(minZ)
            maxZ = autoMax
        if maxZ is not None:
            hist.SetMaximum(maxZ)
        hist.Draw("COLZ,TEXT")
        self.colzPlot = hist
        if textFormat:
            #ROOT.gStyle.SetPaintTextFormat(startingTextFormat)
            pass
        return canv

###############################################################################

def unitTest():
    #create some random data
    r1 = dict()
    colNames = ["col"+str(i) for i in xrange(3)]
    rowNames = ["row"+str(i) for i in xrange(5)]
    rand = ROOT.TRandom3()
    for ir,rn in enumerate(rowNames):
        r1[rn] = [ float(ir+1)+float(ic+1)/10.0 for ic,n in enumerate(colNames) ]
    #create the table
    table = Table.tablefromdict("table_unitTest", r1, headerRow = colNames)
    canv = table.drawColzPlot()
    canv.SaveAs(canv.GetName()+".eps")
    print table.getLatexTable()
    print table.getAsciiTable()
    table.prettyPrint()
    raw_input("wait")
    return

def main():
    unitTest()
    return
    
    
    
    

if __name__ == "__main__":
    main()
    
    
