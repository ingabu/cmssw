from sys import stdout
from ROOT import *

def isEqual(a, b):
    try:
        return a.upper() == b.upper()
    except AttributeError:
        return a == b

def contains(a, b):
    try:
        return b.upper() in a.upper()
    except AttributeError:
        return b in a

##############################################################################

def fillBins(h):
    epsilon=1E-4
    
    if type(h)==type(TH1F()):
        for xBin in range(0,h.GetNbinsX()+2):
            if h.GetBinContent(xBin)<epsilon:
                h.SetBinContent(xBin,epsilon)
                h.SetBinError(xBin,epsilon)
    else:
        for xBin in range(0,h.GetNbinsX()+2):
            for yBin in range(0,h.GetNbinsY()+2):
                if h.GetBinContent(xBin,yBin)<epsilon:
                    h.SetBinContent(xBin,yBin,epsilon)
                    h.SetBinError(xBin,yBin,epsilon)
                
##############################################################################            
            
def unroll(h):
    if type(h)==type(TH1F()): return h

    nX=h.GetNbinsX()
    nY=h.GetNbinsY()

    originalIntegral=h.Integral()

    result=TH1F("", h.GetTitle(), (nX)*(nY), 0, (nX)*(nY))

    for yBin in range(1,nY+1):
        for xBin in range(1,nX+1):
            n = xBin + (nX*(yBin-1))
            result.SetBinContent(n,h.GetBinContent(xBin,yBin))
            result.SetBinError(n,h.GetBinError(xBin,yBin))

    name=h.GetName()
    h.Delete()
    result.SetName(name)

    if abs(originalIntegral-result.Integral())>.001*abs(originalIntegral): raise Exception("utils.py: Unrolled histogram has integral which difffers from original.  Original:",originalIntegral,"Result:",result.Integral())
    return result

##############################################################################

def rebin(h):
    min=0.5
    
    #if type(h)!=type(TH1F()): return h

    #if 

    
##############################################################################
#Printing tables

from math import log10, floor, ceil

def round_sig(x, sig):
    if x==0: return 0

    result=round(x, sig-int(floor(log10(x)))-1)
    if ceil(log10(x)) >= sig: result=int(result)
    return result

def format(number):
    return str(number)
    
def getMaxWidth(table, index):
    #Get the maximum width of the given column index
    max=0
    for row in table:
        try:
            n=len(format(row[index]))
            if n>max: max=n
        except: pass
    return max

def printTable(table,out=stdout):
    """Prints out a table of data, padded for alignment
    @param out: Output stream (file-like object)
    @param table: The table to print. A list of lists.
    Each row must have the same number of columns. """
    col_paddings = []

    maxColumns=0
    for row in table:
        if len(row)>maxColumns: maxColumns=len(row)

    for i in range(maxColumns):
        col_paddings.append(getMaxWidth(table, i))
        
    for row in table:
        # left col
        if row[0]=='break': row[0]='-'*(sum(col_paddings)+(2*len(col_paddings)))
        print >> out, format(row[0]).ljust(col_paddings[0] + 1),
        # rest of the cols
        for i in range(1, len(row)):
            col = format(row[i]).ljust(col_paddings[i] + 2)
            print >> out, col,
        print >> out


##############################################################################

if __name__=='__main__':

    table=[["A","B","C"],[1,2,3],[4,5],[6],['break'],['a long string','short',7,8]]
    printTable(table)
