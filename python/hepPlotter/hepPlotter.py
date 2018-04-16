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
        self.ratio_denominator = False
        self.ratio_numerator   = False
        self.uncertainty       = None
        self.isErrobar = []
        self.yerrors = OrderedDict()




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

        # customizable options
        self.typeOfPlot = typeOfPlot.lower()  # histogram or efficiency plot (ignore capitals)
        self.dimensions = dimensions  # number of dimensions in histogram
        self.ratio_plot = False       # plot a ratio of things
        self.ratio_type = "ratio"     # "ratio","significance"
        self.stacked    = False       # stack plots (1D only)
        self.binning    = 20          # integer for number of bins, or list for non-uniform bins
        self.rebin      = 1           # rebin root histograms
        self.label_size = 20          # size of label for text on plots
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
        self.y_ratio_label  = 'y ratio'
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


    def setValue(self,h,argument,default_value):
        """Set default value if the argument is not set by user"""
        val = None
        try:
            val = kwargs['normed']
        except KeyError: 
            val = default_value

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
        self.setValue( hist,'linestyle','solid' )
        self.setValue( hist,'linewidth',2 )
        self.setValue( hist,'histtype','step' )      # draw_type
        self.setValue( hist,'label','' )
        self.setValue( hist,'yerror',None )
        self.setValue( hist,'weight',None )
        self.setValue( hist,'fmt','o' )
        self.setValue( hist,'mec',hist.linecolor )
        self.setValue( hist,'mfc',hist.color )

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
        """
        if self.dimensions == 1: this_figure = self.plot_hist1d()
        else:                    this_figure = self.plot_hist2d()

        return this_figure



    def plot_hist1d(self):
        """Plot a 1D histogram.  Many options set in 'cms' style file"""
        self.ratio_ylims  = {}
        self.ratio_yticks = {}

        fig      = None
        self.ax1 = None
        self.ax2 = None
        if self.ratio_plot:
            fig = plt.figure()
            gs  = matplotlib.gridspec.GridSpec(2,1,height_ratios=[3,1])
            self.ax1 = fig.add_subplot(gs[0])
            self.ax2 = fig.add_subplot(gs[1],sharex=self.ax1)
            plt.setp(self.ax1.get_xticklabels(),visible=False)

            self.ratio_ylims = {'ymin':{'ratio':0.5,'significance':0.0},
                                'ymax':{'ratio':1.5,'significance':None}}
            self.ratio_yticks = {'ratio':np.asarray([0.6,1.0,1.4]),
                                 'significance':self.ax2.get_yticks()[::2]}
        else:
            fig,self.ax1 = plt.subplots()


        if self.typeOfPlot=="efficiency":
            # draw horizontal lines to guide the eye
            self.ax1.axhline(y=0.25,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=0.50,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=0.75,color='lightgray',ls='--',lw=1,zorder=0)
            self.ax1.axhline(y=1.00,color='lightgray',ls='--',lw=1,zorder=0)

            if self.drawEffDist:
                self.draw_physics_eff()   # Twin axis :: Draw physics distribution with efficiency curve


        ## -- Loop over histograms
        binning     = None
        y_lim_value = None    # weird protection against the axis autoscaling to values smaller than previously drawn histograms
                              # I can't find the source
        max_value   = 0.0
        bottomEdge  = None    # for stacking plots (use this instead of the 'stack' argument
                              # so that all plots can be made in one for-loop
        for n,name in enumerate(self.names):

            h_hist     = self.histograms[name]
            data       = h_hist.data['data']
            error      = h_hist.data['error']
            bin_center = h_hist.data['center']
            bin_width  = h_hist.data['width']
            binning    = h_hist.data['bins']


            if name in self.errorbarplot:
                if self.kwargs[name].get("normed") and self.kwargs[name]["normed"]:
                    # Normalize data to then draw errorbar plot (option not supported by matplotlib)
                    data, bin_edges = np.histogram(bin_center,bins=binning,weights=data,normed=True)

                data  = np.array([i if i else float('NaN') for i in data])  # hide empty values
                NaN_values = np.isnan(data)

                p,c,b = self.ax1.errorbar(bin_center,data,yerr=error,fmt=h_hist.fmt,
                                          mfc=h_hist.mfc,mec=h_hist.mec,label=h_hist.label,
                                          zorder=100,**h_hist.kwarg)
                # record data for later
                data[NaN_values] = 0
                h_hist.plotData  = data
            else:
                if self.draw_types[name]=='step':
                    # Hack for changing legend for step histograms (make a line instead of a rectangle)
                    this_label      = None
                    histStep_pseudo = self.ax1.plot([],[],color=h_hist.linecolor,
                                               lw=h_hist.linewidth,ls=h_hist.linestyle,
                                               label=h_hist.label)
                else:
                    this_label = h_hist.label

                # Make the histogram
                data,b,p = self.ax1.hist(bin_center,bins=binning,weights=data,lw=h_hist.linewidth,
                                 histtype=h_hist.draw_type,bottom=bottomEdge,
                                 ls=h_hist.linestyle,log=h_hist.logplot,color=h_hist.color,
                                 edgecolor=h_hist.linecolor,label=this_label,
                                 **h_hist.kwarg)
                if self.stacked:
                    if bottomEdge is None:
                        bottomEdge  = data
                    else:
                        bottomEdge += data

                # record data for later
                h_hist.plotData = data

                # draw uncertainties -- use GetBinError() or just np.sqrt()
                # only using this for histograms because errorbar has a yerr option
                if self.drawStatUncertainty and self.drawUncertaintyTopFig:
                    self.plotUncertainty(self.ax1,pltname=name,normalize=False)

            # record tallest point for scaling plot
            if max(data) > max_value:
                max_value = max(data)
            if y_lim_value is None or y_lim_value < self.ax1.get_ylim()[1]:   # only increase the autoscale of the y-axis limit
                y_lim_value = self.ax1.get_ylim()[1]
                
        ## End loop over data

        self.binning = np.array(binning) # re-set binning for this instance of hepPlotter

        # y-axis
        if self.ymaxScale is None:
            self.ymaxScale = self.yMaxScaleValues[self.typeOfPlot]
        self.ax1.set_ylim(0.,self.ymaxScale*y_lim_value)
        self.ax1.set_yticks(self.ax1.get_yticks()[1:])
        self.ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        self.setYAxis(self.ax1)

        # x-axis
        self.ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        if self.xlim is not None:
            plt.xlim(self.xlim)

        x_axis = self.ax1
        if self.ratio_plot:
            self.drawRatio()
            x_axis = self.ax2

        self.setXAxis(x_axis)

        # axis ticks
        self.setAxisTickMarks()

        # CMS label
        if self.CMSlabel is not None:
            self.text_labels()

        # Legend
        handles,labels = self.ax1.get_legend_handles_labels() # for re-ordering, if needed
        leg = self.ax1.legend(handles,labels,numpoints=1,
                              ncol=self.numLegendColumns,loc=self.legendLoc)
        leg.draw_frame(False)

        return fig


    def draw_physics_eff(self):
        axTwin = self.ax1.twinx()
        if self.ymaxScale is None:
            self.ymaxScale = self.yMaxScaleValues['efficiency']
        ## Get the histogram
        for effData in self.effData:
            h_effDist = self.hists2plot[effData]
            n,b,p = axTwin.hist(h_effDist['center'],bins=h_effDist['bins'],
                                weights=h_effDist['data'],histtype='step',
                                color=self.colors[effData],lw=self.linewidths[effData],
                                label=self.labels[effData],**self.kwargs[effData])
        axTwin.set_ylabel("",fontsize=0,ha='right',va='top')
        axTwin.set_ylim(ymin=0.0,ymax=self.ymaxScale*max(n))
        # hide y-ticks
        axTwin.set_yticks( np.linspace(axTwin.get_yticks()[0],axTwin.get_yticks()[-1],len(self.ax1.get_yticks())) )
        axTwin.set_yticklabels([])
        self.ax1.set_zorder(axTwin.get_zorder()+1) # put ax in front of axTwin
        self.ax1.patch.set_visible(False)          # hide the 'canvas'

        return


    def setColormap(self):
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

            h_data = np.array(self.hists2plot[key]['data'])
            x_bin_center = self.hists2plot[key]['center']['x']
            y_bin_center = self.hists2plot[key]['center']['y']
            binns_x = self.hists2plot[key]['bins']['x']
            binns_y = self.hists2plot[key]['bins']['y']

                hh,bb = np.histogram2d( x_bin_center,y_bin_center,bins=[binns_x,binns_y],weights=h_data )
                self.colormap = hpt.getDataStructure( hh )

        return


    def plot_hist2d(self):
        """Plot a 2D histogram."""
        fig,self.ax1 = plt.subplots()

        h_data       = np.array(self.hists2plot[key]['data'])
        binns_x      = self.hists2plot[key]['bins']['x']
        binns_y      = self.hists2plot[key]['bins']['y']
        x_bin_center = self.hists2plot[key]['center']['x']
        y_bin_center = self.hists2plot[key]['center']['y']

        # Make the plot -- only one histogram at a time
        # Need to plot ratios in separate histogram
        norm2d = LogNorm() if self.logplot else None
        plt.hist2d(x_bin_center,y_bin_center,bins=[binns_x,binns_y],
                   weights=h_data,cmap=self.colormap,norm=norm2d)

        # Plot bin values, if requested
        if self.bin_yields:
            self.printBinYields()

        # Configure the labels
        self.setYAxis(self.ax1)
        self.setXAxis(self.ax1)
        self.setAxisTickMarks()

        # Configure the colorbar
        cbar = plt.colorbar()
        if self.logplot:
            cbar.ax.set_yticklabels( [r"10$^{\text{%s}}$"%(hpt.extract(i.get_text())) for i in cbar.ax.get_yticklabels()] )

        if self.colorbar_title is not None:
            cbar.ax.set_ylabel(self.colorbar_title)

        # CMS label
        if self.CMSlabel is not None:
            self.text_labels()

        return fig



    def drawRatio(self):
        """
        Ratio plot in frame under main plot
         Can handle plotting multiple ratios with one quantity,
         e.g., compare both up & down systs with nominal
        """
        for i in self.names:
            num_hists       = []   # list of histograms for numerator of ratio
            num_hists_names = []   # names of histograms in numerator of ratio
            if self.ratio_den[i]:
                den_hist = np.array( self.histograms[i] )
                try:
                    num_hists.append( np.array(self.histograms[self.ratio_partner[i]]) )
                    num_hists_names.append( self.ratio_partner[i] )
                except TypeError: # unhashable type: 'list'
                    for j in self.ratio_partner[i]:
                        num_hists.append( np.array(self.histograms[j]) )
                        num_hists_names.append( j )
            else: # only look at the 'ratio_den' terms to prevent plotting ratio twice
                continue


            if self.ratio_type=='ratio':
                self.ax2.axhline(y=1,ls='--',c='k',zorder=1) # draw once


            for nh,num_hist in enumerate(num_hists):

                histName = num_hists_names[nh]

                if self.ratio_type=="ratio":
                    residual     = deepcopy( num_hist / den_hist )
                    residual_err = deepcopy( self.hists2plot[histName]['error'] / den_hist ) # proper uncertainties even for "normalized" plots
                elif self.ratio_type == "significance":
                    residual     = deepcopy( num_hist / np.sqrt(den_hist) )
                    residual_err = [0 for _ in residual] # don't know how to estimate this
                else:
                    print " WARNING :: Un-specified method for ratio plot "
                    print "            Setting ratio equal to 1.0"
                    residual     = np.ones( len(num_hist) )
                    residual_err = [0 for _ in residual]

                # Obtain data points for drawing the ratio
                bin_center = self.hists2plot[histName]['center']
                bin_width  = self.hists2plot[histName]['width']

                if histName in self.errorbarplot:
                    if self.linestyles[histName]=='solid':
                        self.linestyles[histName] = 'o'
                    lcolor = self.linecolors[histName]
                    self.ax2.errorbar(bin_center,residual,xerr=bin_width,yerr=residual_err,
                                  capsize=0,fmt=self.linestyles[histName],
                                  mec=lcolor,mfc=lcolor,color=lcolor,zorder=100)
                else:
                    residual = np.array( [fabs(rr) if (not np.isnan(rr) and not np.isinf(rr)) else 0.0 for rr in residual] )
                    self.ax2.hist(bin_center,bins=self.binning,weights=residual,
                                  edgecolor=self.linecolors[histName],lw=2,
                                  color=self.colors[histName],ls=self.linestyles[histName],
                                  histtype='step',zorder=100)

            ## Simulation Uncertainties
            if self.drawStatUncertainty:
                self.plotUncertainty(self.ax2,pltname=i,normalize=True)

        ## Set the axis properties of the ratio y-axis
        self.ax2.set_ylim(ymin=self.ratio_ylims['ymin'][self.ratio_type],
                          ymax=self.ratio_ylims['ymax'][self.ratio_type])
        self.ratio_yticks['significance']=self.ax2.get_yticks()[::2]

        self.ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        self.ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        self.ax2.set_yticks(self.ratio_yticks[self.ratio_type])
        tickLabels = self.setTickLabels(self.ax2.get_yticks())
        self.ax2.set_yticklabels(tickLabels)
        self.ax2.set_ylabel(self.y_ratio_label,ha='center',va='bottom')

        return



    def plotUncertainty(self,axis,pltname='',normalize=False):
        """
        Plot uncertainties

        @param axis        Which Axes object to draw on
        @param name        Name of sample to plot (access data from global dictionaries)
        @param normalize   draw on ratio plot (divide by total prediction)
        """
        if self.yerrors[pltname] is None:
            return

        error   = self.uncertainties[pltname]
        nominal = self.histograms[pltname]

        # Draw errorbars that are rectangles for each bin
        if self.yerrors[pltname]=='rectangle':
            if normalize:
                resid_unc = {'up': list(((nominal+error)/nominal).repeat(2)),
                             'dn': list(((nominal-error)/nominal).repeat(2))}
            else:
                resid_unc = {'up': list((nominal+error).repeat(2)),
                             'dn': list((nominal-error).repeat(2))}

            fill_between_bins = self.binning   ## for plotting hatch uncertainty
            fill_between_bins = [self.binning[0]]+list(fill_between_bins[1:-1].repeat(2))+[self.binning[-1]]

            axis.fill_between(fill_between_bins,resid_unc['dn'],resid_unc['up'],
                              zorder=10,color=self.colors[pltname],**self.kwargs[pltname])
        # Draw vertical line for errors
        elif self.yerrors[pltname]=='line':
            if normalize:
                resid_unc = {'up': list(((nominal+error)/nominal)),
                             'dn': list(((nominal-error)/nominal))}
            else:
                resid_unc = {'up': list((nominal+error)),
                             'dn': list((nominal-error))}

            error      = [resid_unc['dn'],resid_unc['up']]
            data       = [1. for _ in nominal] if normalize else nominal
            bin_center = 0.5*(self.binning[:-1]+self.binning[1:])
            p,c,b = self.ax1.errorbar(bin_center,data,yerr=error,
                                     fmt=self.linestyles[pltname],color=self.colors[pltname],
                                     zorder=100,**self.kwargs[pltname])

        return


    def setAxis(self,axis,ax=""):
        """Draw labels for a given axis"""
        if not ax: return

        axis_ticks = None
        yaxis = False

        if ax=="y":
            yaxis = True
            axis.set_ylabel(self.y_label,ha='right',va='bottom',position=(0,1))
            axis_ticks = axis.get_yticks()
        elif ax=="x":
            axis.set_xlabel(self.x_label,ha='right',va='top',position=(1,0))
            axis_ticks = axis.get_xticks()
        else: # unsupported option
            return

        # Modify tick labels
        if self.logplot:
            logTickLabels = [r"10$^{\text{%s}}$"%(int(np.log10(i)) ) for i in axis_ticks]
            if yaxis: axis.set_yticklabels(logTickLabels)
            else:     axis.set_xticklabels(logTickLabels)
        else:
            # Draw tick marks as integers or decimal values
            axis_ticks_int = axis_ticks.astype(int)
            if len(set(axis_ticks_int))==len(axis_ticks):
                # all of the ticks are unique as integers so draw them as integers
                if yaxis: axis.set_yticklabels(axis_ticks_int)
                else:     axis.set_xticklabels(axis_ticks_int)
            else:
                # the ticks are not unique as integers, leave them as they are
                axis_ticks = self.setTickLabels(axis_ticks)
                if yaxis: axis.set_yticklabels(axis_ticks[0:-1])
                else:     axis.set_xticklabels(axis_ticks)

        return


    def setYAxis(self,axis):
        """Draw labels for a given axis"""
        self.setAxis(axis,ax="y")
        return


    def setXAxis(self,axis):
        """Draw labels for x-axis"""
        self.setAxis(axis,ax="x")
        return


    def setAxisTickMarks(self):
        """Setup axis ticks"""
        if self.minor_ticks:
            if not self.logplot:
                self.ax1.yaxis.set_minor_locator(self.y1minorLocator) # causes 'tick number error' on logplot
            self.ax1.xaxis.set_minor_locator(self.x1minorLocator)

            if self.ratio_plot:
                self.ax2.xaxis.set_minor_locator(self.x2minorLocator)
                self.ax2.yaxis.set_minor_locator(self.y2minorLocator)

        return



    def text_labels(self):
        """
        Write the CMS, Energy, and LUMI labels
        """
        if self.dimensions==2 and self.CMSlabel!='outer':
            print " WARNING :: You have chosen a label position "
            print "            not considered for 2D plots. "
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



    def printBinYields(self):
        """Print bin values inside the plots for each bin."""
        h_data       = None
        x_bin_center = None
        y_bin_center = None
        binns_x      = None
        binns_y      = None

        for key in self.names:
            if key in self.effData: continue

            h_data = np.array(self.hists2plot[key]['data'])
            x_bin_center = self.hists2plot[key]['center']['x']
            y_bin_center = self.hists2plot[key]['center']['y']
            binns_x = self.hists2plot[key]['bins']['x']
            binns_y = self.hists2plot[key]['bins']['y']

        N = len(x_bin_center) # number of bins (a flattened list for both TH1/TH2)

        text_colors = self.bin_yields_color
        if text_colors is None:
            text_colors = ['k' for _ in range(N)]

        for i in range(N):
            self.ax1.text(x_bin_center[i],y_bin_center[i],
                          "{0:.1f}".format(h_data[i]),
                          ha='center',va='center',
                          color=text_colors[i])

        return


    def setTickLabels(self,values):
        """Remove extra trailing characters"""
        values = [str(i) for i in values]
        len_ticklabels = min( [len(i) for i in values if len(i)>0] )
        values = [i[:len_ticklabels] for i in values]
        return values


    def savefig(self):
        """Save the figure"""
        plt.savefig(self.saveAs+'.'+self.format,format=self.format)
        plt.close()

        return



## THE END
