"""
Created:        6 April     2016
Last Updated:  15 April     2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

Bennett Magy
bmagy@umichSPAMNOT.edu
University of Michigan, Ann Arbor, MI 48109
-----

Class to make a simple instance each time we want some basic plots!

This does not include an interface to load/access data.
Here we just plot the data we're given.

Base class for turning histograms or efficiency curves into plots
"""
import os
import sys
import ROOT
from math import fabs
from copy import deepcopy
from collections import OrderedDict

import numpy as np
import matplotlib as mpl
mpl_version = mpl.__version__
mpl.style.use('cms')

from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LogNorm
from matplotlib.ticker import AutoMinorLocator,FormatStrFormatter

import hepPlotterTools as hpt
import hepPlotterLabels as hpl



class Histogram(object):
    """Class for containing histogram objects to plot"""
    def __init__(self):
        self.name  = ''
        self.color = 'k' 
        self.linecolor = 'k'
        self.linestyle = '-'
        self.linewidth = 2
        self.draw_type = 'hist'
        self.label  = ''
        self.weight = None
        self.data   = None
        self.stack  = False              # ability to stack some plots, not all
        self.ratio_denominator = False
        self.ratio_numerator   = False
        self.uncertainty       = None
        self.isErrobar  = False  # plt.errorbar()
        self.isLinePlot = False  # basic plt.plot()
        self.yerrors    = None
        self.plotData   = None   # the data prepared for the plot (maybe slightly modified)



class HepPlotter(object):
    def __init__(self,typeOfPlot,dimensions):
        """
        @param typeOfPlot    Set the kind of plot: histogram or efficiency
        @param dimensions    Number of dimension for the histogram/efficiency: 1 or 2
        """
        if not isinstance(dimensions,(int,long)):
            print " You have specified an unsupported type for 'dimension'"
            print " For the hepPlotter class, choose either 1 or 2 dimenions."
            print " Exiting. "
            sys.exit(1)

        # customizable options for the figure
        self.typeOfPlot = typeOfPlot.lower()  # histogram or efficiency plot (ignore capitals)
        self.dimensions = dimensions  # number of dimensions in histogram
        self.ratio_plot = False       # plot a ratio of things
        self.ratio_type = "ratio"     # "ratio","significance"
        self.stacked    = False       # stack plots (1D only)
        self.binning    = 20          # integer for number of bins, or list for non-uniform bins
        self.rebin      = 1           # rebin root histograms
        self.normed     = False       # normalize histogram
        self.logplot    = False       # plot on log scale
        self.underflow  = False       # plot the underflow
        self.overflow   = False       # plot the overflow
        self.colormap   = None        # 2D plot colormap
        self.colorbar_title = None
        self.xlim       = None        # tuple for (xmin,xmax)
        self.ymaxScale  = None        # scale y-axis
        self.format     = 'pdf'          # file format for saving image
        self.saveAs     = "plot_{0}d_{1}".format(self.dimensions,self.CMSlabelStatus) # save figure with name
        self.bin_yields = False       # print bin yields inside histogram
        self.bin_yields_color = None  # array of text colors for each bin
        self.drawEffDist    = False    # draw the physics distribution for efficiency (jet_pt for jet trigger)
        self.x_label        = 'x'
        self.y_label        = 'y'
        self.y_ratio_label  = 'ratio'
        self.extra_text     = hpl.PlotText()
        self.minor_ticks    = True
        self.lumi           = '14.7'
        self.plotLUMI       = False
        self.CMSlabel       = 'top left'     # 'top left', 'top right' & 'outer' for 2D
        self.CMSlabelStatus = 'Internal'  # ('Simulation')+'Internal' || 'Preliminary'
        self.numLegendColumns = 2
        self.legendLoc        = 1
        self.drawStatUncertainty   = False
        self.drawUncertaintyTopFig = False  # draw uncertainties in the top frame
        self.uncertaintyHistType   = 'step'

        self.text_coords = {'top left': {'x':[0.03]*3,'y':[0.97,0.90,0.83]},\
                            'top right':{'x':[0.97]*3,'y':[0.97,0.90,0.83]},\
                            'outer':    {'x':[0.02,0.99,0.99],'y':[1.0,1.0,0.9]}}

        # Arguments for plotting uncertainties [color = kGreen-8 , alpha=0.5]
        self.p_hatch_args    = {'hatch':'','color':'#99cc99','edgecolor':'none','alpha':0.5}
        self.yMaxScaleValues = {"histogram":1.3,"efficiency":1.3}


        self.ratio_ylims = {'ymin':{'ratio':0.5,'significance':0.0},
                            'ymax':{'ratio':1.5,'significance':None}}
        self.ratio_yticks = {'ratio':np.asarray([0.6,1.0,1.4]),
                             'significance':[0,1,2]}

        return


    def initialize(self):
        """Initialize some things."""
        self.ax1          = None
        self.ax2          = None
        self.kwargs       = {}

        self.histograms   = OrderedDict()    # {'name',Histogram()}

        if self.format!='pdf': print " WARNING : Chosen format '{0}' may conflict with PDF backend"

        ## 2D plot
        if self.dimensions==2: self.setColormap()

        # draw minor ticks in the 'right' places
        self.x1minorLocator = AutoMinorLocator()
        self.y1minorLocator = AutoMinorLocator()
        self.x2minorLocator = AutoMinorLocator()
        self.y2minorLocator = AutoMinorLocator()
        self.yTwinMinorLocator = AutoMinorLocator()   # twin axis for efficiency plots

        return


    def setValue(self,h,argument,default_value,arg_other=''):
        """Set default value if the argument is not set by user"""
        val = None
        try:
            val = kwargs[argument]
        except KeyError: 
            if arg_other:
                try: 
                    val = kwargs[arg_other]
                except KeyError: 
                    val = default_value
            else: val = default_value

        setattr(h,argument,val)

        return


    def setDefaults(self,hist,**kwargs):
        """Set some default values for the plots if the aren't set by the user
           This also ensures matplotlib defaults are used where necessary
        """
        user_kwargs = kwargs.keys()
        for kw in user_kwargs:
            setattr( hist,kw,kwargs[kw] )       # set user-defined coords

        # Set default attributes (if not set by the user)
        self.setValue( hist,'normed',False )
        self.setValue( hist,'linecolor','k' )
        self.setValue( hist,'color','k' )
        self.setValue( hist,'linestyle','solid','ls' )
        self.setValue( hist,'linewidth',2,'lw' )
        self.setValue( hist,'draw_type','step' )
        self.setValue( hist,'label','' )
        self.setValue( hist,'yerror',None )
        self.setValue( hist,'weight',None )
        self.setValue( hist,'fmt','o' )
        self.setValue( hist,'mec',hist.linecolor,'markeredgecolor' )
        self.setValue( hist,'mfc',hist.color,'markerfacecolor' )
        self.setValue( hist,'isErrorbar',False )
        self.setValue( hist,'isLinePlot',False )
        self.setValue( hist,'stack',self.stacked )

        return


    def Add(self,data,name='',weights=None,ratio_den=False,ratio_num=False,ratio_partner=None,**kwargs):
        """
        Add histogram data for this figure.

        @param data             data for plot (python array or ROOT TH1)
        @param name             name to identify histogram object
        @param weights          weights for making histogram data
        @param ratio_den        denominator in ratio plot - True/False
        @param ratio_num        numerator of ratio plot   - True/False
        @param ratio_partner    Object to take ratio of with this one
        @param kwargs           arguments for matplotlib objects (histogram 1/2D; errorbar)
                   -- hist:     https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist.html
                   -- hist2d:   https://matplotlib.org/api/_as_gen/matplotlib.pyplot.hist2d.html
                   -- errorbar: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.errorbar.html
                   -- lineplot: https://matplotlib.org/api/_as_gen/matplotlib.pyplot.plot.html
        """
        hist = Histogram()
        hist.name = name

        # Internally process data -- convert to relevant format for plotting (ROOT->matplotlib)
        if self.typeOfPlot=="histogram" and isinstance(data,ROOT.TH1):
            # Make histogram plot with TH1/TH2
            if self.dimensions==1:
                h_data = hpt.hist2list(data,name=name,reBin=self.rebin,normed=hist.normed)
            else:
                h_data = hpt.hist2list2D(data,name=name,reBin=self.rebin,normed=hist.normed)
        elif self.typeOfPlot=="efficiency" and isinstance(data,ROOT.TEfficiency):
            # Make efficiency plot with TEfficiency
            h_data = hpt.TEfficiency2list(data)
        elif self.typeOfPlot=="efficiency" and isinstance(data,ROOT.TH1):
            # Make efficiency plot -- draw a histogram
            h_data = hpt.hist2list(data,name=name,reBin=self.rebin,normed=True)
            self.drawEffDist = True
        else:
            # Catch-all (Numpy data)
            if self.dimensions==1:
                h_data = hpt.data2list(data,weights=weights,normed=hist.normed,binning=self.binning)
            else:
                h_data = hpt.data2list2D(data,weights=weights,normed=hist.normed,binning=self.binning)

        setDefaults(hist,kwargs)

        hist.ratio_den  = ratio_den          # denominator in ratio
        hist.ratio_num  = ratio_num          # numerator in ratio
        hist.ratio_partner = ratio_partner   # other object to plot against in ratio
        hist.data = h_data

        # store in this list to use throughout the class
        self.histograms[name] = hist

        return


    def execute(self):
        """
        Execute the plot.
        return the Figure object to the user (they can edit it if they please)
        Many options set in 'cms' style file
        """
        # Setup the figure
        fig      = None
        self.ax1 = None
        self.ax2 = None

        if self.ratio_plot:
            fig = plt.figure()
            gs  = matplotlib.gridspec.GridSpec(2,1,height_ratios=[3,1])
            self.ax1 = fig.add_subplot(gs[0])
            self.ax2 = fig.add_subplot(gs[1],sharex=self.ax1)
            plt.setp(self.ax1.get_xticklabels(),visible=False)
        else:
            fig,self.ax1 = plt.subplots()


        ## -- Efficiency specific items
        if self.typeOfPlot=="efficiency":
            # draw horizontal lines to guide the eye
            self.ax1.axhline(y=0.25,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=0.50,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=0.75,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=1.00,color='lightgray',ls='--',lw=1,zorder=0)

            if self.drawEffDist:
                self.draw_physics_eff()   # Twin axis :: Draw physics distribution with efficiency curve


        ## -- Loop over objects to plot
        binning     = None
        y_lim_value = None    # protection against the axis autoscaling to values smaller than previously drawn histograms
        max_value   = 0.0
        bottomEdge  = None    # for stacking plots (use this instead of the 'stack' argument
                              # so that all plots can be made in one for-loop)

        for n,name in enumerate(self.names):
            h_hist     = self.histograms[name]

            if h_hist.isErrorbar:
                self.draw_errorbar_plot(h_hist)
            elif h_hist.isLinePlot:
                self.draw_line_plot(h_hist)
            else:
                if self.dimensions==1:
                    self.draw_histogram(h_hist)
                else:
                    self.draw_histogram2D(h_hist)

            # record tallest point for scaling plot
            this_max = max(h_hist.data['data'])
            if this_max > max_value:
                max_value = this_max
            if y_lim_value is None or y_lim_value < self.ax1.get_ylim()[1]:   # only increase the autoscale of the y-axis limit
                y_lim_value = self.ax1.get_ylim()[1]


        ## AXES
        # y-axis
        if self.ymaxScale is None:
            self.ymaxScale = self.yMaxScaleValues[self.typeOfPlot]
        self.ax1.set_ylim(ymax=self.ymaxScale*y_lim_value)
        self.ax1.set_yticks(self.ax1.get_yticks()[1:])
        #self.ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.setAxis(self.ax1,ax="y",label=self.y_label)

        # x-axis
        #self.ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if self.xlim is not None:
            plt.xlim(self.xlim)

        x_axis = self.ax1
        if self.ratio_plot:
            self.drawRatio()
            x_axis = self.ax2

        self.setAxis(x_axis,ax="x",label=self.x_label)

        # axis ticks
        self.setMinorAxisTickMarks()

        ## CMS && Text labels
        if self.CMSlabel is not None or len(self.extra_text.texts)>0:
            self.writeText()

        ## Legend
        self.drawLegend()

        return fig


    def drawLegend(self):
        """Draw legend"""
        handles,labels = self.ax1.get_legend_handles_labels() # for re-ordering, if needed
        leg = self.ax1.legend(handles,labels,ncol=self.numLegendColumns,loc=self.legendLoc)  # no ax2 legend, for now

        return


    def getAxes(self):
        """Return the axes to the user"""
        if self.ratio_plot:
            return self.ax1,self.ax2
        else:
            return self.ax1


    def draw_line_plot(self,h):
        """Draw basic plt.plot()"""
        self.ax1.plot(h.data['center'],h.data['data'],**h.kwarg)

        return


    def draw_errorbar_plot(self,h):
        """Draw errorbar plot"""
        if h.normed:
            # Normalize data using hist then draw errorbar plot (option not supported by matplotlib)
            data, bin_edges = np.histogram(bin_center,bins=binning,weights=data,normed=True)

        data = np.array([i if i else float('NaN') for i in data])  # hide empty values
        NaN_values = np.isnan(data)

        p,c,b = self.ax1.errorbar(bin_center,data,yerr=error,fmt=h_hist.fmt,
                                  mfc=h_hist.mfc,mec=h_hist.mec,label=h_hist.label,
                                  zorder=100,**h_hist.kwarg)
        # record data for later
        data[NaN_values] = 0
        h_hist.plotData  = data

        return


    def draw_histogram(self,h,axis=None):
        """Draw histogram"""
        if h.draw_type=='step':
            # Hack for changing legend for step histograms (make a line instead of a rectangle)
            this_label      = None
            histStep_pseudo = self.ax1.plot([],[],color=h_hist.linecolor,lw=h_hist.linewidth,ls=h_hist.linestyle,label=h_hist.label)
        else:
            this_label = h_hist.label

        # Make the histogram
        data,b,p = axis.hist(bin_center,bins=binning,weights=data,lw=h_hist.linewidth,
                             histtype=h_hist.draw_type,bottom=bottomEdge,
                             ls=h_hist.linestyle,log=h_hist.logplot,color=h_hist.color,
                             edgecolor=h_hist.linecolor,label=this_label,
                             **h_hist.kwargs)
        if self.stacked or h.stack:
            # raise the bottom edge for the next plot to stack on this one
            if bottomEdge is None:
                bottomEdge  = data
            else:
                bottomEdge += data

        # record data for later
        h.plotData = data

        # draw uncertainties -- use GetBinError()
        # only using this for histograms because errorbar has a yerr option
        if self.drawStatUncertainty and self.drawUncertaintyTopFig:
            self.plotUncertainty(self.ax1,pltname=name,normalize=False)

        return


    def draw_physics_eff(self):
        """Draw the corresponding physics distribution with an efficiency curve using twin axis"""
        axTwin = self.ax1.twinx()
        if self.ymaxScale is None:
            self.ymaxScale = self.yMaxScaleValues['efficiency']
        # Get the histogram
        for effData in self.effData:
            h_effDist = self.histograms[effData]
self.draw_histogram(h_effDist,axTwin)

            n,b,p = axTwin.hist(h_effDist['center'],bins=h_effDist['bins'],
                                weights=h_effDist['data'],histtype='step',
                                color=self.colors[effData],lw=self.linewidths[effData],
                                label=self.labels[effData],**self.kwargs[effData])
        # Suppress labels
kwargs = {"fontsize":0}
setAxis(axTwin,ax="y",kwargs)

        axTwin.set_ylabel("",fontsize=0,ha='right',va='top')
        axTwin.set_ylim(ymin=0.0,ymax=self.ymaxScale*max(n))
        # hide y-ticks
        axTwin.set_yticks( np.linspace(axTwin.get_yticks()[0],axTwin.get_yticks()[-1],len(self.ax1.get_yticks())) )
        axTwin.set_yticklabels([])

        self.ax1.set_zorder(axTwin.get_zorder()+1) # put ax in front of axTwin
        self.ax1.patch.set_visible(False)          # hide the 'canvas'

        return


    def setColormap(self,h):
        """Colormap setup for 2D plots"""
        linear_cmap_choice  = np.random.choice(["Reds","Blues","Greens"])
        default_cmap_choice = np.random.choice(["viridis","magma","inferno","plasma"])

        colormaps = {'diverge':"bwr",
                     'linear':linear_cmap_choice,
                     'default':default_cmap_choice}
        try:
            self.colormap = getattr( plt.cm,colormaps[self.colormap] )   # use a pre-defined choice
        except AttributeError:
            try:
                self.colormap = getattr(plt.cm,self.colormap)            # access map from matplotlib choices
            except AttributeError:
                print " WARNING : Unsupported colormap '{0}'".format(self.colormap)
                print "           Choosing the colormap based on data structure "

                h_data       = h.data['data']
                x_bin_center = h.data['center']['x']
                y_bin_center = h.data['center']['y']
                binns_x      = h.data['bins']['x']
                binns_y      = h.data['bins']['y']

                hh,bb = np.histogram2d( x_bin_center,y_bin_center,bins=[binns_x,binns_y],weights=h_data )
                self.colormap = hpt.getDataStructure( hh )

        return


    def draw_histogram2D(self,h,axis=None):
        """Plot a 2D histogram."""
        h_data       = h.data['data']
        binns_x      = h.data['bins']['x']
        binns_y      = h.data['bins']['y']
        x_bin_center = h.data['center']['x']
        y_bin_center = h.data['center']['y']

        # Make the plot -- only one histogram at a time
        # Need to plot ratios in separate histogram
        norm2d = LogNorm() if self.logplot else None
        plt.hist2d(x_bin_center,y_bin_center,bins=[binns_x,binns_y],
                   weights=h_data,cmap=self.colormap,norm=norm2d)

        # Plot bin values, if requested
        if self.bin_yields:
            self.printBinYields()   # only supported for 2D plots

        # Configure the colorbar
        cbar = plt.colorbar()
     setAxis(cbar,ax='y')

        if self.logplot: cbar.ax.set_yticklabels( [r"10$^{\text{%s}}$"%(hpt.extract(i.get_text())) for i in cbar.ax.get_yticklabels()] )
        if self.colorbar_title is not None: cbar.ax.set_ylabel(self.colorbar_title)

        return



    def drawRatio(self):
        """
        Ratio plot in frame under main plot
         Can handle plotting multiple ratios with one quantity,
         e.g., compare both up & down systs with nominal
        """
        if self.ratio_type=='ratio':
            self.ax2.axhline(y=1,ls='--',c='k',zorder=1) # draw once

        for i in self.histograms:
            num_hists       = []   # list of histograms for numerator of ratio
            num_hists_names = []   # names of histograms in numerator of ratio

            if i.ratio_den:
                # This data is the denominator in ratio_plot
                den_hist = i.plotData
                try:
                    num_hists.append( np.array(self.histograms[i.ratio_partner]) )
                    num_hists_names.append( i.ratio_partner )
                except TypeError: # unhashable type: 'list'
                    for j in i.ratio_partner:
                        num_hists.append( np.array(self.histograms[j]) )
                        num_hists_names.append( j )
            else: # only look at the 'ratio_den' terms to prevent plotting ratio twice
                continue


            for nh,num_hist in enumerate(num_hists):

                histName = num_hists_names[nh]

                if self.ratio_type=="ratio":
                    residual     = deepcopy( num_hist.plotData / den_hist )
                    residual_err = deepcopy( num_hist.data['error'] / den_hist ) # proper uncertainties even for "normalized" plots
                elif self.ratio_type == "significance":
                    residual     = deepcopy( num_hist.plotData / np.sqrt(den_hist) )
                    residual_err = [0 for _ in residual] # don't know how to estimate this
                else:
                    print " WARNING :: Un-specified method for ratio plot '{0}'. Setting ratio equal to 1.0".format(self.ratio_type)
                    residual     = np.ones( len(num_hist.data['center']) )
                    residual_err = [0 for _ in residual]

                # Obtain data points for drawing the ratio
                bin_center = num_hist.data['center']
                bin_width  = num_hist.data['width']

                if num_hist.isErrorbar:

draw errorbar plot
                    self.ax2.errorbar(bin_center,residual,xerr=bin_width,yerr=residual_err,
                                      fmt=num_hist.linestyle,mec=num_hist.mec,mfc=num_hist.mfc,zorder=100,**num_hist.kwarg)

draw line plot

                else:

draw histogram
                    residual = np.array( [fabs(rr) if (not np.isnan(rr) and not np.isinf(rr)) else 0.0 for rr in residual] )
                    self.ax2.hist(bin_center,bins=self.binning,weights=residual,
                                  edgecolor=num_hist.linecolor,lw=num_hist.linewidth,
                                  color=num_hist.color,ls=num_hist.linestyle,
                                  histtype=num_hist.draw_type,zorder=100,**num_hist.kwarg)

            ## Simulation Uncertainties
            if self.drawStatUncertainty:
                self.plotUncertainty(self.ax2,pltname=i.name,normalize=True)

        ## Set the axis properties of the ratio y-axis
        self.ax2.set_ylim(ymin=self.ratio_ylims['ymin'][self.ratio_type],
                          ymax=self.ratio_ylims['ymax'][self.ratio_type])
        self.ratio_yticks['significance']=self.ax2.get_yticks()[::2]
        #self.ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        #self.ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        self.ax2.set_yticks(self.ratio_yticks[self.ratio_type])
        self.setAxis(self.ax2,ax="y",label=self.y_ratio_label)

        return


    def plotUncertainty(self,axis,pltname='',normalize=False):
        """
        Plot uncertainties

        @param axis        Which Axes object to draw on
        @param name        Name of sample to plot (access data from global dictionaries)
        @param normalize   draw on ratio plot (divide by total prediction)
        """
        if not pltname: return

        h_hist  = self.histograms[pltname]
        error   = h_hist.data['error']
        nominal = h_hist.plotData
        binning = h_hist.data['binning']

        # Draw errorbars that are rectangles for each bin
        if h_hist.yerror=='rectangle':
            if normalize:
                resid_unc = {'up': list(((nominal+error)/nominal).repeat(2)),
                             'dn': list(((nominal-error)/nominal).repeat(2))}
            else:
                resid_unc = {'up': list((nominal+error).repeat(2)),
                             'dn': list((nominal-error).repeat(2))}

            fill_between_bins = binning   ## for plotting hatch uncertainty
            fill_between_bins = [binning[0]]+list(fill_between_bins[1:-1].repeat(2))+[binning[-1]]

            axis.fill_between(fill_between_bins,resid_unc['dn'],resid_unc['up'],
                              zorder=10,color=h_hist.color,**h_hist.kwarg)
        # Draw vertical line for errors
        elif h_hist.yerror=='line':
            if normalize:
                resid_unc = {'up': list(((nominal+error)/nominal)),
                             'dn': list(((nominal-error)/nominal))}
            else:
                resid_unc = {'up': list((nominal+error)),
                             'dn': list((nominal-error))}

            error      = [resid_unc['dn'],resid_unc['up']]
            data       = [1. for _ in nominal] if normalize else nominal
            bin_center = 0.5*(binning[:-1]+binning[1:])

draw errorbar plot

            p,c,b = self.ax1.errorbar(bin_center,data,yerr=error,
                                     fmt=h_hist.fmt,color=h_hist.color,
                                     zorder=100,**h_hist.kwarg)

        return


    def setAxis(self,axis,ax="",label='',ratio=False):
        """Draw labels for a given x-/y-axis"""
        if not ax: return
        if ax!="y" and ax!="x": return

        axis_ticks = None
        if ax=="y":
            axis.set_ylabel(label,ha='right',va='bottom',position=(0,1))
            axis_ticks = axis.get_yticks()
        else:
            axis.set_xlabel(label,ha='right',va='top',position=(1,0))
            axis_ticks = axis.get_xticks()

        # Modify tick labels and write in scientific notation (only for the upper frame)
        if self.logplot and not ratio:
            logTickLabels = [r"10$^{\text{%s}}$"%(int(np.log10(i)) ) for i in axis_ticks]
            if ax=='y': axis.set_yticklabels(logTickLabels)
            else:       axis.set_xticklabels(logTickLabels)

        return


    def setMinorAxisTickMarks(self):
        """Setup axis ticks"""
        if not self.minor_ticks: return

        if not self.logplot:
            self.ax1.yaxis.set_minor_locator(self.y1minorLocator) # causes 'tick number error' on logplot
        self.ax1.xaxis.set_minor_locator(self.x1minorLocator)

        if self.ratio_plot:
            self.ax2.xaxis.set_minor_locator(self.x2minorLocator)
            self.ax2.yaxis.set_minor_locator(self.y2minorLocator)

        return


    def writeText(self):
        """Write the CMS, Energy, & LUMI labels + any extra text the user wants
           This only writes text on self.ax1, nothing on self.ax2
        """
        if self.dimensions==2 and self.CMSlabel!='outer':
            print " WARNING :: You have chosen a label position "
            print "            unsupported for 2D plots. Unknown consequences! "
            print "            Please consider changing the "
            print "            parameter 'CMSlabel' to 'outer'."

        text = self.text_coords[self.CMSlabel]

        cms_stamp    = hpl.CMSStamp(self.CMSlabelStatus)
        lumi_stamp   = hpl.LumiStamp(self.lumi)
        energy_stamp = hpl.EnergyStamp()

        cms_stamp.coords    = [text['x'][0], text['y'][0]]
        lumi_stamp.coords   = [text['x'][1], text['y'][1]]  # not used right now, always drawn with the energy
        energy_stamp.coords = [text['x'][1], text['y'][1]]

        # modify defaults
        if self.CMSlabel == 'top right':
            cms_stamp.ha    = 'right'    # move text labels appropriately
            lumi_stamp.ha   = 'right'
            energy_stamp.ha = 'right'
        if self.dimensions==2:
            cms_stamp.va    = 'bottom'   # change alignment for 2d labels
            lumi_stamp.ha   = 'right'
            energy_stamp.ha = 'right'


        self.ax1.text(cms_stamp.coords[0],cms_stamp.coords[1],cms_stamp.text,fontsize=cms_stamp.fontsize,
                      ha=cms_stamp.ha,va=cms_stamp.va,transform=self.ax1.transAxes)

        energy_lumi_text = energy_stamp.text+", "+lumi_stamp.text if self.plotLUMI else energy_stamp.text
        self.ax1.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_lumi_text,
                      fontsize=energy_stamp.fontsize,ha=energy_stamp.ha,va=energy_stamp.va,
                      color=energy_stamp.color,transform=self.ax1.transAxes)


        ## Extra text -- other text labels the user wants to add
        for txtItem in self.extra_text.texts:
            if txtItem.transform is None: txtItem.transform = self.ax1.transAxes
            self.ax1.text( txtItem.coords[0],txtItem.coords[1],txtItem.text,
                           fontsize=txtItem.fontsize,ha=txtItem.ha,va=txtItem.va,
                           color=txtItem.color,transform=txtItem.transform  )

        return



    def printBinYields(self,h):
        """Print bin values inside the plots for each bin."""
        h_data       = None
        x_bin_center = None
        y_bin_center = None
        binns_x      = None
        binns_y      = None

        for key in self.names:
            h_data = h.data['data']
            x_bin_center = h.data['center']['x']
            y_bin_center = h.data['center']['y']
            binns_x = h.data['bins']['x']
            binns_y = h.data['bins']['y']

        text_colors = self.bin_yields_color
        if text_colors is None: text_colors = ['k' for _ in x_bin_center]  # default option

        for i,bc in enumerate(x_bin_center):
            self.ax1.text(bc,y_bin_center[i],"{0:.1f}".format(h_data[i]),
                          ha='center',va='center',color=text_colors[i])

        return


    def savefig(self):
        """Save the figure"""
        plt.savefig(self.saveAs+'.'+self.format,format=self.format)
        plt.close()

        return



## THE END
