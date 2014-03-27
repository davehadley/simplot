import ROOT
import StringIO
import pprint
import itertools
import os

#import prettytable

###############################################################################

class Table:
    def __init__(self, name):
        self.name = name
        self.setHeaderRow(None)
        self.rows = []
        self.setFloatingPointFormatString("{:.2g}")
        
    def setFloatingPointFormatString(self, formatStr):
        self.floatingPointFormatString = formatStr
    
    def setHeaderRow(self, headerRow):
        self.headerRow = None
        #ensure that header row is all strings
        if headerRow is not None:
            self.headerRow = [str(e) for e in headerRow]
        
    def addRow(self, row):
        self.rows.append(row)
    
    def getNRows(self):
        nRows = len(self.rows)
        return nRows
    
    def getNCols(self):
        nCols = max( (len(l) for l in self.rows) )
        return nCols
    
    def hasRowLabels(self):
        return any((len(r)>0 and type(r[0]) is str for r in self.rows))
        
    @staticmethod
    def tableFromDictionary(name, container, headerRow=None):
        table = Table(name)
        if headerRow:
            table.setHeaderRow(headerRow)
        for key in sorted(container.keys()):
            row = container[key]
            inputRow = [key] + row
            table.addRow(inputRow)
        return table
    
    @staticmethod
    def tableFromMatrix(name, matrix, xLabels=None, yLabels=None):
        nRows = matrix.GetNrows()
        nCols = matrix.GetNcols()
        if xLabels is None:
            xLabels = [str(i) for i in range(nCols)]
        if yLabels is None:
            yLabels = [str(i) for i in range(nRows)]
        tab = Table(name)
        tab.setHeaderRow(yLabels)
        for j in xrange(nRows):
            row = [ xLabels[j] ]
            for i in xrange(nCols):
                v = matrix[i][j]
                row.append(v)
            tab.addRow(row)
        return tab
    
    def writeLatexTable(self, fileName):
        try:
            os.makedirs(os.path.dirname(fileName))
        except os.error:
            pass
        outFile = open(fileName,"w")
        print >>outFile,self.getLatexTable()
        return 
    
    def getLatexTable(self):
        latex = StringIO.StringIO()
        nRows = self.getNRows()
        nCols = self.getNCols()
        newLine = " \\\\\n"
        colFormat = "|".join(["c"]*nCols)
        print >>latex,"{"
        print >>latex,"%\\tiny" 
        print >>latex,"\\begin{tabular}{"+colFormat+"}"
        allData = [self._getPrintedHeaderRow()]+self.rows
        for row in allData:
            strRow = self._convertRowToString(row)
            strRow = self._sanitiseLatexRow(strRow)
            line = " & ".join(strRow) + newLine
            print >>latex,line,
        print >>latex,"\\end{tabular}"
        print >>latex,"}" 
        return latex.getvalue()
    
    def getAsciiTable(self):
        pt = prettytable.PrettyTable(self._getPrintedHeaderRow())
        for row in self.rows:
            row = [ self._convertEntryToString(n) for n in row]
            pt.add_row(row)
        return str(pt)
    
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
    
    def prettyString(self):
        stream = StringIO.StringIO()
        print >>stream,"---- Table(",self.name,") ----"
        headerRow = self._getPrintedHeaderRow()
        allData = [headerRow]
        for row in self.rows:
            allData.append( self._convertRowToString(row) )
        pprint.pprint(allData, stream=stream, indent=2)
        return stream.getvalue()
    
    def prettyPrint(self):
        print self.prettyString()
        
    def __str__(self):
        return "Table("+self.name+")"
    
    def _getPrintedHeaderRow(self):
        result = None
        if len(self.headerRow)==self.getNCols()-1:
            result = [ "" ] + self.headerRow
        else:
            result = list(self.headerRow)
        return result
    
    def _convertRowToString(self, row):
        strRow = [self._convertEntryToString(v) for v in row]
        return strRow

    def _convertEntryToString(self, entry):
        strValue = None
        if isinstance(entry, basestring):
            #already a string, copy value
            strValue = str(entry)
        else:
            #not a string, try to get value and format
            try:
                value,error = self._getValueAndError(entry)
                vs = self.floatingPointFormatString.format(value)
                es = ""
                if error is not None:
                    es = self.floatingPointFormatString.format(error)
                    strValue = vs + " +- " + es
                else:
                    strValue = vs
            except:
                #can't convert to string, fall back on standard python string conversion
                strValue = str(entry)
        return strValue
    
    def _getValueAndError(self, entry):
        value,error = None,None
        try:
            #is it a single number?
            value = float(entry)
        except:
            try:
                #try unpacking a tuple
                value,error = entry
                value = float(value)
                error = float(error)
            except:
                #don't know what to do now, raise an exception
                raise Exception("can't convert entry to value and error",entry)
        return value,error
    
    def _sanitiseLatexRow(self, row):
        result = [self._sanitiseLatexString(entry) for entry in row]
        return result
    
    def _sanitiseLatexString(self, entry):
        result = entry
        #deal with underscores
        result = result.replace("\\_","_")
        result = result.replace("_","\\_")
        #change \pm
        result = result.replace("+-","\\ensuremath{\\pm}")
        return result
    


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
    table = Table.tableFromDictionary("table_unitTest", r1, headerRow = colNames)
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
    
    
