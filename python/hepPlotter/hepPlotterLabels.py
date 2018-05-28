"""
Created:        --
Last Updated:   16 February  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

File for containing information about plotting.
Note: in hepPlotter package, function "hist1d" is deprecated
"""
from array import array
from collections import OrderedDict
from hepPlotterTools import hist1d


class Sample(object):
    """Class for organizing plotting information about physics samples"""
    def __init__(self,label='',color=''):
        self.label = label
        self.color = color

class Variable(object):
    """Class for organizing plotting information about variables"""
    def __init__(self,binning=[],label=''):
        self.binning = binning
        self.label   = label


def sample_labels():
    """Dictionaries that contain Samples objects."""
    ## Sample information
    samples = {}

    # Standard Model
    ttbar = r't$\bar{\text{t}}$'
    samples['ttbar']     = Sample(label=ttbar,color='white')
    samples['dijet']     = Sample(label=r'Dijets', color='purple')
    samples['multijet']  = Sample(label=r'Multi-jet', color='purple')
    samples['diboson']   = Sample(label=r'Diboson',color='green')
    samples['singletop'] = Sample(label=r'Single Top',color='blue')
    samples['ttbarW']    = Sample(label=ttbar+'W',color='#C9FFE5')
    samples['ttbarZ']    = Sample(label=ttbar+'Z',color='#7FFFD4')
    samples['ttbarV']    = Sample(label=ttbar+'V',color='cyan')
    samples['ttbarH']    = Sample(label=ttbar+'H',color='#3AB09E')
    samples['ttbarX']    = Sample(label=ttbar+'V',color='#008B8B')
    samples['wjets']     = Sample(label=r'W+jets',color='yellow')
    samples['zjets']     = Sample(label=r'Z+jets',color='darkorange')

    samples['ttbar_QONLY'] = Sample(label=ttbar+" (q)",    color='#696969')
    samples['ttbar_BONLY'] = Sample(label=ttbar+" (b)",    color='#808080')
    samples['ttbar_BQ']    = Sample(label=ttbar+" (qb)",   color='#79160F')
    samples['ttbar_W']     = Sample(label=ttbar+" (qq)",   color='#DC682A')
    samples['ttbar_FULL']  = Sample(label=ttbar+" (qqb)",  color='red')
    samples['ttbar_OTHER'] = Sample(label=ttbar+" (other)",color='#EC706B')
    samples['ttbar_NONE']  = Sample(label=ttbar+" (other)",color='#696969')

    # Machine learning (CWoLa)
    samples['top'] = Sample(label='Top',color='white')
    samples['antitop'] = Sample(label='Anti-top',color='white')

    # Machine Learning (AK8+AK4)
    samples['QB'] = Sample(label=ttbar+' (QB)',color='white')
    samples['W']  = Sample(label=ttbar+' (W)',color='white')

    # Data
    samples['data']      = Sample(label=r'Data',color='black')

    # Signal
    samples['signal']      = Sample(label='Signal',color='Reds')
    samples['zprime_1000'] = Sample(label=r'm$_{\text{Z}^\prime}$=1.0 TeV',color='r')

    # Generic
    samples['mu']        = Sample(label=r'$\mu$+jets',color='k')
    samples['mujets']    = Sample(label=r'$\mu$+jets',color='k')
    samples['el']        = Sample(label=r'e+jets',color='k')
    samples['ejets']     = Sample(label=r'e+jets',color='k')
    samples['muel']      = Sample(label=r'$\ell$+jets',color='k')
    samples['ljets']     = Sample(label=r'$\ell$+jets',color='k')

    return samples


def variable_labels():
    """Dictionaries that contain Variables objects."""
    _phi  = r'$\phi$'
    _eta  = r'$\eta$'
    _T    = r'$_\text{T}$ [GeV]'
    _mass = 'Mass [GeV]'

    variables = {}

    variables['ljet_C2']    = Variable(binning=hist1d(10,  0.,   0.6), label=r'Large-R Jet C$_2^{\beta\text{=1}}$')
    variables['ljet_D2']  = Variable(binning=hist1d(20,  0.,   5.0), label=r'Large-R Jet D$_2^{\beta\text{=1}}$')
    variables['ljet_d12'] = Variable(binning=hist1d(20,  0.,  125.), label=r'Large-R Jet $\sqrt{\text{d}_{\text{12}}}$ [GeV]')
    variables['ljet_d23'] = Variable(binning=hist1d(12,  0.,   60.), label=r'Large-R Jet $\sqrt{\text{d}_{\text{23}}}$ [GeV]')
    variables['ljet_eta'] = Variable(binning=hist1d(20, -3.,    3.), label=r'Large-R Jet '+_eta)
    variables['ljet_phi'] = Variable(binning=hist1d(20, -2.,    2.), label=r'Large-R Jet $\phi$')
    variables['ljet_m']   = Variable(binning=hist1d(40,  0.,  400.), label=r'Large-R Jet '+_mass)
    variables['ljet_pt']  = Variable(binning=hist1d(14,200., 1500.), label=r'Large-R Jet p'+_T)
    variables['ljet_BEST_t'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(top)')
    variables['ljet_BEST_w'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(W)')
    variables['ljet_BEST_z'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(Z)')
    variables['ljet_BEST_h'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(H)')
    variables['ljet_BEST_j'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(jet)')
    variables['ljet_BEST_t_j'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet BEST(top/(top+jet))')
    variables['ljet_tau1'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet $\tau_{\text{1}}$')
    variables['ljet_tau2'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet $\tau_{\text{2}}$')
    variables['ljet_tau3'] = Variable(binning=hist1d(10,  0.,   1.0), label=r'Large-R Jet $\tau_{\text{3}}$')
    variables['ljet_tau21'] = Variable(binning=hist1d(11,  0.,   1.1), label=r'Large-R Jet $\tau_{\text{21}}$')
    variables['ljet_tau32'] = Variable(binning=hist1d(11,  0.,   1.1), label=r'Large-R Jet $\tau_{\text{32}}$')
    variables['ljet_charge'] = Variable(binning=hist1d(80,  -5.,   5.),label=r'Large-R Jet Charge')
    variables['ljet_SDmass'] = Variable(binning=hist1d(40,0.,400.),    label=r'Large-R Jet Softdrop Mass [GeV]')
    variables['ljet_softDropMass'] = Variable(binning=hist1d(40,0.,400.), label=r'Large-R Jet Softdrop Mass [GeV]')

    for i in range(16):
        variables['ljet_deepAK8_{0}'.format(i)] = Variable(binning=hist1d(10,0,1), label=r'Large-R Jet DeepAK8[{0}]'.format(i))
    variables['ljet_jet_m']      = Variable(binning=hist1d(50,0.,5000.), label=r'Large-R Jet + Small-R Jet '+_mass)
    variables['ljet_jet_deltaR'] = Variable(binning=hist1d(10,0.,5.),    label=r'$\Delta$R(Large-R Jet,Small-R Jet)')

    variables['ljet_subjet_0_charge_Qpos'] = Variable(binning=hist1d(50,-5,5), label=r'Large-R Jet Subjet 0 charge')
    variables['ljet_subjet_0_bdisc_Qpos']  = Variable(binning=hist1d(10, 0,1), label=r'Large-R Jet Subjet 0 b-disc.')
    variables['ljet_subjet_1_charge_Qpos'] = Variable(binning=hist1d(50,-5,5), label=r'Large-R Jet Subjet 1 charge')
    variables['ljet_subjet_1_bdisc_Qpos']  = Variable(binning=hist1d(10, 0,1), label=r'Large-R Jet Subjet 1 b-disc.')
    variables['ljet_subjet_0_charge_Qneg'] = Variable(binning=hist1d(50,-5,5), label=r'Large-R Jet Subjet 0 charge')
    variables['ljet_subjet_0_bdisc_Qneg']  = Variable(binning=hist1d(10, 0,1), label=r'Large-R Jet Subjet 0 b-disc.')
    variables['ljet_subjet_1_charge_Qneg'] = Variable(binning=hist1d(50,-5,5), label=r'Large-R Jet Subjet 1 charge')
    variables['ljet_subjet_1_bdisc_Qneg']  = Variable(binning=hist1d(10, 0,1), label=r'Large-R Jet Subjet 1 b-disc.')

    variables['ljet_subjet0_charge'] = Variable(binning=hist1d(50,-5,5),  label=r'Large-R Jet Subjet 0 charge')
    variables['ljet_subjet0_bdisc']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 0 b-disc.')
    variables['ljet_subjet0_mass']   = Variable(binning=hist1d(20, 0,200),label=r'Large-R Jet Subjet 0 '+_mass)
    variables['ljet_subjet0_mrel']   = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 0 Relative '+_mass)
    variables['ljet_subjet0_ptrel']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 0 Relative p'+_T)
    variables['ljet_subjet0_tau21']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 0 $\tau_{\text{21}}$')
    variables['ljet_subjet0_tau32']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 0 $\tau_{\text{32}}$')
    variables['ljet_subjet1_charge'] = Variable(binning=hist1d(50,-5,5),  label=r'Large-R Jet Subjet 1 charge')
    variables['ljet_subjet1_bdisc']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 1 b-disc.')
    variables['ljet_subjet1_mass']   = Variable(binning=hist1d(20, 0,200),label=r'Large-R Jet Subjet 1 '+_mass)
    variables['ljet_subjet1_mrel']   = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 1 Relative '+_mass)
    variables['ljet_subjet1_ptrel']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 1 Relative p'+_T)
    variables['ljet_subjet1_tau21']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 1 $\tau_{\text{21}}$')
    variables['ljet_subjet1_tau32']  = Variable(binning=hist1d(10, 0,1),  label=r'Large-R Jet Subjet 1 $\tau_{\text{32}}$')
    variables['ljet_subjets_deltaQ'] = Variable(binning=hist1d(100,-10,10),  label=r'$\Delta$Q (Large-R Jet Subjets)')
    variables['ljet_contain']        = Variable(binning=hist1d(11,-5.5,5.5), label=r'Large-R Jet Containment')

    variables['jet_pt']  =   Variable(binning=hist1d(10,  25., 500), label=r'Small-R Jet p'+_T)
    variables['jet_eta'] =   Variable(binning=hist1d(10, -2.5, 2.5), label=r'Small-R Jet '+_eta)
    variables['btags_n'] =   Variable(binning=hist1d( 4, -0.5, 3.5), label=r'Number of b-tags')
    variables['n_btags'] =   Variable(binning=hist1d( 6, -0.5, 5.5), label=r'Number of b-tags')
    variables['n_jets']  =   Variable(binning=hist1d(11, -0.5,10.5), label=r'Number of Small-R Jets')
    variables['n_ljets'] =   Variable(binning=hist1d( 6, -0.5, 5.5), label=r'Number of Large-R Jets')
    variables['jet_bdisc'] = Variable(binning=hist1d(10,   0.,  1.), label=r'Small-R Jet b-disc.')

    variables['lep_eta'] = Variable(binning=hist1d(10,-2.5,   2.5),label=r'Lepton '+_eta)
    variables['lep_pt']  = Variable(binning=hist1d(10, 25.,  300.),label=r'Lepton p'+_T)

    variables['mu_pt']  = Variable(binning=hist1d(10, 25.,  300.),label=r'Muon p'+_T)
    variables['mu_eta'] = Variable(binning=hist1d(10,-2.5,   2.5),label=r'Muon '+_eta)
    variables['mu_phi'] = Variable(binning=hist1d(10,-2.5,   2.5), label=r'Muon '+_phi)
    variables['mu_ptrel'] = Variable(binning=hist1d(10,0,100), label=r'Muon p$_\text{T}^\text{rel}$')
    variables['mu_drmin'] = Variable(binning=hist1d(10,0,  5), label=r'Muon $\Delta$R$_\text{min}$')
    variables['el_pt']  = Variable(binning=hist1d(10, 25.,  300.),label=r'Electron p'+_T)
    variables['el_eta'] = Variable(binning=hist1d(10,-2.5,   2.5),label=r'Electron '+_eta)
    variables['el_phi'] = Variable(binning=hist1d(10,-2.5,   2.5), label=r'Electron '+_phi)
    variables['el_ptrel'] = Variable(binning=hist1d(10,0,100), label=r'Electron p$_\text{T}^\text{rel}$')
    variables['el_drmin'] = Variable(binning=hist1d(10,0,  5), label=r'Electron $\Delta$R$_\text{min}$')

    variables['deltaR_lep_ak4'] = Variable(binning=hist1d(10,0,  5), label=r'$\Delta$R(lepton,AK4)')
    variables['pTrel_lep_ak4']  = Variable(binning=hist1d(10,0,100), label=r'p$_\text{T}^\text{rel}$(lepton,AK4)')
    variables['deltaR_lep_ak8'] = Variable(binning=hist1d(10,0,  5), label=r'$\Delta$R(lepton,AK8)')
    variables['deltaR_ak4_ak8'] = Variable(binning=hist1d(10,0,  5), label=r'$\Delta$R(AK4,AK8)')

    variables['nu_phi']  = Variable(binning=hist1d(64,-3.2,  3.2), label=r'$\nu$ '+_phi)
    variables['nu_eta']  = Variable(binning=hist1d(10,-2.5,  2.5), label=r'$\nu$ '+_eta)
    variables['nu_pt']   = Variable(binning=hist1d(30,   0,  600.),  label=r'$\nu$ p'+_T)

    variables['ht']      = Variable(binning=hist1d(50,  0., 5000.), label=r'H'+_T)
    variables['HT']      = Variable(binning=hist1d(50,  0., 5000.), label=r'H'+_T)
    variables['mtw']     = Variable(binning=hist1d(12,  0.,  120.), label=r'm$_\text{T}^\text{W}$ [GeV]')
    variables['mlb']     = Variable(binning=hist1d(32,  0.,  800.), label=r'm$_{\ell\text{b}}$')
    variables['mass_lb'] = Variable(binning=hist1d(32,  0.,  800.), label=r'm$_{\ell\text{b}}$')
    variables['met_met'] = Variable(binning=hist1d(50,  0.,  500.), label=r'E$_{\text{T}}^{\text{miss}}$ [GeV]')
    variables['met_phi'] = Variable(binning=hist1d( 6,-3.2,   3.2), label=r'$\phi^{\text{miss}}$ [GeV]')

    ttbar = r"\text{t}\bar{\text{t}}"
    variables['mtt']  = Variable(binning=hist1d(25,0,5000),label=r'm$_{%s}$ [GeV]'%ttbar)
    variables['pttt'] = Variable(binning=hist1d(10,0,500),label=r'p$_{\text{T},%s }$ [GeV]'%ttbar)
    variables['ytt']  = Variable(binning=hist1d(10,0,5),label=r'y$_{%s}$ [GeV]'%ttbar)
    variables['beta'] = Variable(binning=hist1d(10,0,1),label=r'$\beta_{z,%s}$ [GeV]'%ttbar)
    variables['dy']   = Variable(binning=hist1d(12,-3,3),label=r'$\Delta|\text{y}|$')

    variables['deltay']    = variables['dy']
    variables['mttbar']    = variables['mtt']
    variables['pTttbar']   = variables['pttt']
    variables['yttbar']    = variables['ytt']
    variables['betatt']    = variables['beta']
    variables['betattbar'] = variables['beta']

    return variables




## -- Classes for handling text on plots
class Text(object):
    """Class to hold extra text object"""
    def __init__(self):
        self.text     = ''         # Actual text to show on plot
        self.coords   = [0.,0.]    # coordinates on plot to draw to the text
        self.fontsize = 16
        self.color    = 'k'
        self.ha = 'left'           # horizontal alignment
        self.va = 'top'            # vertical alignment
        self.transform = None      # 'None' so the user can change it -- it will be set below
        return

    def __str__(self):
        """print text object with attributes"""
        for i in ['text','coords','fontsize','color','ha','va','transform']:
            print "%-*s: %s" % (10,i,self.__dict__[i])
        return


class PlotText(object):
    """Class to draw new text on the plots"""
    def __init__(self):
        self.texts  = []

    def Add(self,plt_text,**txt_kwargs):
        """
	Add new text to the plot
        @param plt_text    text to draw
        @param txt_kwargs  key-word arguments:
                             'coords','fontsize','color','ha','va','transform'
        """
	pltTextObject = Text()
        pltTextObject.text = plt_text

        # set parameters of the text object if they are passed to the 'Add()' function
        # - use defaults if no argument is passed - this ensures any unsupported 
        # arguments don't harm anything in the text object

        for param in dir(pltTextObject):
            if param.startswith("__"): continue
            try:
                setattr( pltTextObject,param,txt_kwargs[param] )
            except KeyError: # use the defaults
                continue

        self.texts.append(pltTextObject)

        return

    def Print(self):
        """Print out the text arguments"""
        for text in self.texts:
            print text
    def getText(self):
        """Return the list of Text objects"""
        return self.texts

class EnergyStamp(Text):
    """Class for writing center of mass energy on plot"""
    def __init__(self):
        Text.__init__(self)
        self.text = r"(13 TeV)"
        self.fontsize = 18
        self.ha = 'right'
        self.va = 'bottom'

class LumiStamp(Text):
    """Class for writing luminosity on plot"""
    def __init__(self,lumi="36.1"):
        Text.__init__(self)
        self.text = r"%s fb$^{\text{-1}}$"%(lumi)
        self.fontsize = 18
        self.ha = 'right'
        self.va = 'bottom'

class CMSStamp(Text):
    """Class for writing official CMS name & plot type (Simulation, Internal, etc.) on plot"""
    def __init__(self,label_status="Internal"):
        Text.__init__(self)
        self.text = r"\textbf{CMS} {\Large \textit{%s}}"%(label_status)
        self.fontsize = 18
        self.ha = 'left'
        self.va = 'top'



if __name__ == '__main__':
    print "Do not execute this file, only import it."


## THE END ##
