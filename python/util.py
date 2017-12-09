"""
Created:         1 August 2017
Last Updated:    5 August 2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

File that holds any and all misc. functions 
to be called from other python scripts.
(All information in one file => one location to update)
"""
import ROOT



class LorentzVector(object):
    """Simple class that extends TLorentzVectors"""
    def __init__(self):
        self.p4 = ROOT.TLorentzVector()
        self.charge = 0
        self.mv2c10 = 0
        self.mv2c10_bin = 0
        self.index = -1
        return


class Ljet(LorentzVector):
    """Simple class for containing all top quark information"""
    def __init__(self):
        LorentzVector.__init__(self)

        self.target = -1       # target value for NN training
        self.tjets  = []       # list of track jets ghost-associated to ljet
        self.n_tjets = 0       # Number of ghost-associated track jets
        self.reco_mttbar  = 0  # Reconstructed m_ttbar
        self.truth_mttbar = 0  # Generator-level m_ttbar
        return

    def __str__(self):
        """Print properties for debugging"""
        command  = " > Top Quark Jet "
        command += "\n   - Target            = {0}".format(self.target)
        command += "\n   - Index             = {0}".format(self.index)
        command += "\n   - Charge            = {0}".format(self.charge)

        return command




def deltaR(tlvA,tlvB,dR=0.75):
    """Determine if two objects are matched in eta-phi space"""
    result = False
    try:
        result = (tlvA.p4.DeltaR(tlvB.p4)<dR)
    except:
        return False # something that isn't a LorentzVector

    return result



def file2list(filename):
    """Load text file and dump contents into a list"""
    listOfFiles = open( filename,'r').readlines()
    listOfFiles = [i.rstrip('\n') for i in listOfFiles]

    return listOfFiles



class VERBOSE(object):
    """Object for handling output"""
    def __init__(self):
        self.verboseMap = {"DEBUG":0,
                           "INFO": 1,
                           "WARNING":2,
                           "ERROR":  3};
        self.level = "WARNING"

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

    def verbose(self,level,message):
        if self.verboseMap[level] >= self.verboseMap[self.level]:
            print " {0} :: {1}".format(level,message)
        return


## THEN ##
