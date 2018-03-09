"""
Created:        --
Last Updated:    2 March 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

File that holds any and all misc. functions 
to be called from other python scripts.
(All information in one file => one location to update!)
"""
import ROOT
import numpy as np


def getHistSeparation( S, B ):
    """Compare TH1* S and B -- need same dimensions
       Copied from : https://root.cern.ch/doc/master/MethodBase_8cxx_source.html#l02740
    """
    separation = 0

    nstep  = S.GetNbinsX()
    xaxis  = S.GetXaxis()

    nS = S.GetSumOfWeights()
    nB = B.GetSumOfWeights()
    for bin in range(nstep):
        s = S.GetBinContent( bin+1 )/nS
        b = B.GetBinContent( bin+1 )/nB
        if (s+b)>0: separation += (s - b)*(s - b)/(s + b)

    separation *= 0.5

    return separation


def GetSeparation2D( S, B ):
    """Compare TH2* S and B -- need same dimensions"""
    separation = 0

    nbinsx = S.GetNbinsX()
    xaxis  = S.GetXaxis()

    nbinsy = S.GetNbinsY()
    yaxis  = S.GetYaxis()

    integral_s = S.Integral()
    integral_b = B.Integral()

    for x in range(nbinsx):
        for y in range(nbinsy):
            s = S.GetBinContent( x+1,y+1 )/integral_s
            b = B.GetBinContent( x+1,y+1 )/integral_b

            if (s+b) > 0: separation += (s - b)*(s - b)/(s + b)

    separation *= 0.5

    return separation



def getSeparation(sig,bkg):
    """Calculate separation between two distributions"""
    separation = 0

    nS = np.sum(sig)
    nB = np.sum(bkg)
    for ss,bb in zip(sig,bkg):
        s = ss/nS
        b = bb/nB
        
        if (s+b) > 0: separation += (s - b)*(s - b)/(s + b)

    separation *= 0.5

    return separation


def read_config(filename,separation=" "):
    """
    Read configuration file with data stored like:
       'config option'
    And the 'config' and 'option' are separated by a character, e.g., " "
    """
    data = file2list(filename)
    cfg = {}
    for i in data:
        j = i.split(separation)
        cfg[j[0]] = j[1]
    return cfg


def extract(str_value, start_='{', stop_='}'):
    """Extract a string between two symbols, e.g., parentheses."""
    extraction = str_value[str_value.index(start_)+1:str_value.index(stop_)]
    return extraction


def file2list(filename):
    """Load text file and dump contents into a list"""
    listOfFiles = open( filename,'r').readlines()
    listOfFiles = [i.rstrip('\n') for i in listOfFiles if not i.startswith("#")]
    return listOfFiles


def str2bool(param):
    """Convert a string to a boolean"""
    if param in ['true','True','1']:
        return True
    else:
        return False


class VERBOSE(object):
    """Object for handling output"""
    def __init__(self):
        self.verboseMap = {"DEBUG":0,
                           "INFO": 1,
                           "WARNING":2,
                           "ERROR":  3};
        self.level     = "WARNING"
        self.level_int = 2

    def initialize(self):
        """Setup the integer level value"""
        self.level_int = self.verboseMap[self.level]

    def level_value(self):
        """Return the integer value"""
        return self.level_int

    def DEBUG(self,message):
        """Debug level - most verbose"""
        self.verbose("DEBUG",message)
        return

    def INFO(self,message):
        """Info level - standard output"""
        self.verbose("INFO",message)
        return

    def WARNING(self,message):
        """Warning level - if something seems wrong but code can continue"""
        self.verbose("WARNING",message)
        return

    def ERROR(self,message):
        """Error level - something is wrong"""
        self.verbose("ERROR",message)
        return

    def compare(self,level1,level2=None):
        """Compare two levels"""
        if level2 is None:
            return self.verboseMap[level1]>=self.level_int
        else:
            return self.verboseMap[level1]>=self.verboseMap[level2]

    def verbose(self,level,message):
        """Print message to the screen"""
        if self.compare( level ):
            print " {0} :: {1}".format(level,message)
        return

## THE END ##