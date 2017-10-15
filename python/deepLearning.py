"""
Created:        11 November  2016
Last Updated:   11 December  2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----

Class for performing deep learning to classify 
large-R jets as either a top or anti-top quark.

https://keras.io/
"""
import json
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import hepPlotterLabels as hpl
import hepPlotterTools as hpt
from hepPlotter import HepPlotter
from keras import __version__ as vs
from keras.models import Sequential,model_from_json
from keras.layers import Dense, Activation
from keras.utils.visualize_util import plot
from sklearn.cross_validation import train_test_split,StratifiedKFold
from sklearn.metrics import roc_curve, auc

from matplotlib import rc
rc('font', family='sans-serif')



# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

class DeepLearning(object):
    """Deep Learning for Top/Anti-Top discrimination"""
    def __init__(self):
        # event-level features included -- target '1' = leading top; '0' = leading antitop
        self.features = ['charge',
                         'subjets_btag_charge',\
                         'subjets_nbtag_charge',\
                         'subjets_delta_btag_nbtag']#,\
                         #'subjet_pt0_charge',\
                         #'subjet_pt0_btag_charge',\
                         #'subjet_pt0_nbtag_charge']
        self.features2plot = self.features+['pt','eta',#'reco_m_ttbar',\
                                            'subjet_pt0_btag_charge',\
                                            'subjet_pt0_nbtag_charge']
                              # 'qb1qnb2_qb2qnb1','deltaQ'
        self.text_dicts = {'pt':    {'label':r'Large-R Jet p$_\text{T}$ [GeV]',\
                                     'bins':hpl.hist1d(40,200,2000)},\
                           'eta':   {'label':r'Large-R Jet $\eta$',\
                                     'bins':hpl.hist1d(12,-3,3)},\
                           'charge':{'label':r'Large-R Jet Charge [e]',\
                                     'bins':hpl.hist1d(96,-2,2)},\
                           'charge_other':{'label':r'Large-R Jet Charge [e]',\
                                     'bins':hpl.hist1d(96,-3,3)},\
                           'subjets_btag_charge': {'label':r'Track Subjet Charge [e] (b-tagged {0}\%)',\
                                                   'bins':hpl.hist1d(24,-3,3)},\
                           'subjets_nbtag_charge':{'label':r'Track Subjet Charge [e] (non b-tagged {0}\%)',\
                                                   'bins':hpl.hist1d(24,-3,3)},\
                           'subjet_pt0_charge':   {'label':r'Leading Track Subjet Charge [e]',\
                                                   'bins':hpl.hist1d(24,-3,3)},\
                           'subjet_pt0_btag_charge':   {'label':r'Leading B-tagged Track Subjet Charge [e]',\
                                                   'bins':hpl.hist1d(24,-3,3)},\
                           'subjet_pt0_nbtag_charge':   {'label':r'Leading Not B-tagged Track Subjet Charge [e]',\
                                                   'bins':hpl.hist1d(24,-3,3)},\
                           'reco_m_ttbar':{'label':r'm$_{\text{t}\bar{\text{t}}}$',\
                                           'bins':hpl.hist1d(47,350,5050)},\
                           'data':{'label':'Test','bins':[ (i-4) for i in range(20)] },\
                           'qb1qnb2_qb2qnb1':{'label':r'(Q$_\text{b1}$+Q$_\text{!b2}$)-(Q$_\text{b2}$+Q$_\text{!b1}$)',
                                              'bins':hpl.hist1d(96,-3,3)},\
                           'deltaQ':{'label':'Top Charge - Antitop Charge [e]','bins':hpl.hist1d(96,-3,3)},\
                           'subjets_delta_btag_nbtag':{'label':r'$\Delta$Q(b-tag,non-b-tag)','bins':hpl.hist1d(96,-3,3)}}
        self.processlabel_args = {'text':r"G$_\text{RS}\rightarrow\text{t}\bar{\text{t}}$",
                                  'coords':[0.03,0.83]}
                                  # r"Z$^\prime$(m=2.25TeV)$\rightarrow$ t$\bar{\text{t}}$"

        self.history    = None
        self.date       = datetime.date.today().strftime('%d%b%Y')
        self.hep_data   = None
        self.dnn_data   = None
        self.score      = None # set later
        self.df         = None # set later
        self.model      = None # set later
        self.fpr        = None # set later
        self.tpr        = None # set later
        self.output     = 'data/DNN/' # set later
        self.metadata   = {}   # set later
        self.train_data = {}   # set later
        self.test_data  = {}   # set later
        self.dnn_name   = "ttbar_DNN"  # name to access in lwtnn ('variables.json')

        ## NN parameters
        self.batch_size = 32
        self.nb_epoch   = 1        
        self.loss       = 'binary_crossentropy' # preferred for binary classification
        self.optimizer  = 'sgd'
        self.metrics    = ['accuracy']
        self.init       = 'normal' #'uniform'
        self.input_dim  = len(self.features)
        self.output_dim = 3*len(self.features)
        self.nHiddenLayers = 1

        return



    def getHEPData(self):
        """
        Load the data for NN.  
        This uses data created by root2keras.py where there is metadata 
        saved to the json output (to reproduce results later &
        understand how they were produced) <- remove these extra keys.
        """
        data  = json.load( open(self.hep_data,'r') )
        fdata = dict(  (key,value) for key,value in data.items() if key in self.features2plot+['target'] )
        self.metadata = data['metadata']

        # -- convert to DataFrame for easier slicing 
        self.df = pd.DataFrame( fdata )
        self.df.hist()

        # -- split the dataset into train and test portions
        X_data  = self.df[self.features].values
        Y_data  = self.df['target'].values
        X_train, X_test, y_train, y_test = train_test_split(X_data,Y_data,test_size=0.4)

        self.test_data  = {'X':X_test, 'Y':y_test}
        self.train_data = {'X':X_train,'Y':y_train}

        return



    def buildNN(self):
        """Initialize the NN"""
        self.model = Sequential() # The Keras Sequential model is a linear stack of layers.

        if not self.train_data and not self.test_data:
            self.get_data()

        # Add 1st layer
        self.model.add( Dense(self.output_dim, input_dim=self.input_dim, init=self.init) )
        self.model.add( Activation("relu") )

        # Add middle layer(s) [identical setup for now]
        # Number of middle layers = hidden layers - 1 (initial layer above)
        for _ in range(self.nHiddenLayers-1):
            self.model.add( Dense(self.output_dim, init=self.init) )
            self.model.add( Activation("relu") )

        # Add the output
        self.model.add( Dense(1,init=self.init) )
        self.model.add( Activation("sigmoid") )

        # Build the model
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)

        # define 10-fold cross validation test harness
        cvscores = []
        X = self.train_data['X']
        Y = self.train_data['Y']
        kfold = StratifiedKFold(Y, n_folds=10, shuffle=True, random_state=seed)
        for ind,(train,test) in enumerate(kfold):
            print "   - Fitting K-Fold {0}".format(ind)
            self.model.fit(X[train], Y[train],nb_epoch=self.nb_epoch, 
                           batch_size=self.batch_size, verbose=0)
            # evaluate the model
            print "     + Evaluate the model: "
            scores = self.model.evaluate(X[test], Y[test], verbose=0,
                                         batch_size=self.batch_size)
            print "       {0}: {1:.2f}%".format(self.model.metrics_names[1], scores[1]*100)
            cvscores.append(scores[1] * 100)

        print "   Finished K-Fold cross-validation: "
        print "   - Accuracy: {0:.2f}% (+/- {1:.2f}%)".format(np.mean(cvscores), np.std(cvscores))


#        self.history = self.model.fit(self.train_data['X'],self.train_data['Y'],
#                                      nb_epoch=self.nb_epoch,batch_size=self.batch_size)

        return



    def saveModel(self):
        """Save the model for use later"""
        output = self.output+'/'+self.hep_data.split('/')[-1].split('.')[0]+'_'+self.date

        ## model architecture
        json_string = self.model.to_json()
        with open(output+'_model.json', 'w') as outfile:
            outfile.write(json_string)

        ## save the weights of a model, you can do so in HDF5
        self.model.save_weights(output+'_weights.h5')

        ## Plot the model to view it later
        plot(self.model,to_file=output+'_model.png',show_shapes=True)

        return



    def loadModel(self):
        """Load existing model to make plots or predictions"""
        json_file  = open(self.dnn_data+"_model.json", 'r')
        model_json = json_file.read()

        self.model = model_from_json(model_json)
        self.model.load_weights(self.dnn_data+"_weights.h5")

        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)

        return



    def getScore(self):
        """Return the score from a test sample"""
        self.score = self.model.predict( self.test_data['X'], batch_size=self.batch_size )

        return self.score


    def plot_features(self):
        """Plot the features"""
        print " Plotting features comparing top quarks and anti-quarks. "
        top  = self.df.loc[self.df['target'] == self.metadata['t_target']]
        tbar = self.df.loc[self.df['target'] == self.metadata['tbar_target']]

        processed_features = []
        for hi,feature in enumerate(self.features2plot):

            eventlevel = False
            if feature.startswith('t_'):      # specific top properties in dataframe
                feature = feature[2:]
                eventlevel = True
            elif feature.startswith('tbar_'): # specific tbar properties in dataframe
                feature = feature[5:]
                eventlevel = True
            else:
                eventlevel = False # single object in dataframe (use 'target' to distinguish)

            if feature in processed_features: continue
            else: processed_features.append(feature)

            if 'btag' in feature:
                x_label = self.text_dicts[feature]['label'].format( self.metadata['btag_wkpt'] )
            else:
                x_label = self.text_dicts[feature]['label']

            hist = HepPlotter("histogram",1)

            hist.ratio_plot  = False
            hist.binning     = self.text_dicts[feature]['bins']
            hist.stacked     = False
            hist.logplot     = False
            hist.x_label     = x_label
            hist.y_label     = "Events"
            hist.format      = 'png'
            hist.saveAs      = self.output+"/hist_"+feature+"_"+self.date
            hist.ATLASlabel       = 'top left'
            hist.ATLASlabelStatus = 'Simulation Internal'
            hist.numLegendColumns = 1
            hist.extra_text.Add(self.processlabel_args['text'],**self.processlabel_args)

            hist.initialize()

            multiply=1.
            if feature=='reco_m_ttbar' or feature=='pt':
                multiply = 1e-3

            if eventlevel:
                hist.Add(self.df['t_'+feature],name=feature+'_top',linecolor='r',color='r',
                         draw='step',label='Large-R Jet (top)')
                hist.Add(self.df['tbar_'+feature],name=feature+'_tbar',linecolor='b',color='b',
                         draw='step',label='Large-R Jet (anti-top)')
            else:
                hist.Add(top[feature].multiply(multiply),name=feature+'_top',linecolor='r',color='r',
                         draw='step',label='Large-R Jet (top)')
                hist.Add(tbar[feature].multiply(multiply),name=feature+'_tbar',linecolor='b',color='b',
                         draw='step',label='Large-R Jet (anti-top)')

            p = hist.execute()
            hist.savefig()


        ## Correlation Matrices of Features (top/antitop) ##
        corrmat_df_top  = top[self.features].corr()
        corrmat_df_tbar = tbar[self.features].corr()

        names = ["top","tbar"]
        fontProperties = {'family':'sans-serif'}
        opts = {'cmap': plt.get_cmap("bwr"), 'vmin': -1, 'vmax': +1}

        for c,corrmat in enumerate([corrmat_df_top,corrmat_df_tbar]):
            fig,ax = plt.subplots()

            heatmap1 = ax.pcolor(corrmat, **opts)
            cbar = plt.colorbar(heatmap1, ax=ax)

            cbar.ax.set_yticklabels( [i.get_text().strip('$') for i in cbar.ax.get_yticklabels()], **fontProperties )

            labels = corrmat.columns.values
            labels = [i.replace('_','\_') for i in labels]
            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_xticklabels(labels, fontProperties, fontsize=18, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(labels, fontProperties, fontsize=18, minor=False)

            text_args = {'fontsize':18,'ha':'left','va':'bottom','transform':ax.transAxes}

            ## ATLAS Label + Signal name
            ax.text(0.02,1.00,r"\textbf{\textit{ATLAS}} Simulation Internal",**text_args)
            ax.text(0.03,0.93, r"G$_\text{RS}\rightarrow\text{t}\bar{\text{t}}$",**text_args)

            ## Energy Label
            text_args['ha'] = 'right'
            ax.text(0.99,1.00,r"$\sqrt{\text{s}}$ = 13 TeV",**text_args)

            plt.savefig(self.output+"/correlations_"+names[c]+"_"+self.date+".png",
                        format='png',dpi=300,bbox_inches='tight')
            plt.close()

        return


    def save_features(self):
        """
        Save the features to a json file to load in the lwtnn later
        Hard-coded scale & offset; must change later if necessary
        """
        text = """  {
    "inputs": ["""
        for feature in self.features:
            tmp = """
      {"name": "%(feature)s",
       "scale":  1,
       "offset": 0},""" % {'feature':feature}
            text += tmp
        text += "],"
        text += """
    "class_labels": ["%(name)s"],
    "keras_version": "%(version)s",
    "miscellaneous": {}
  }
""" % {'version':vs,'name':self.dnn_name}
        varsFileName = self.output+'/variables.json'
        varsFile     = open(varsFileName,'w')
        varsFile.write(text)

        return


    def plot_score(self):
        """Plot the features"""

        train_scores = self.model.predict( self.train_data['X'], batch_size=self.batch_size )
        #train_scores = self.model.evaluate( self.train_data['X'],self.train_data['Y'] )
        top_train_scores  = train_scores[ self.train_data['Y']==1 ]
        tbar_train_scores = train_scores[ self.train_data['Y']==0 ]
        top_test_scores   = self.score[ self.test_data['Y']==1 ]
        tbar_test_scores  = self.score[ self.test_data['Y']==0 ]

        hist = HepPlotter("histogram",1)

        hist.ratio_plot  = False
        hist.normed      = True
        hist.binning     = [0.05*i for i in range(21)]
        hist.stacked     = False
        hist.logplot     = False
        hist.x_label     = "DNN Score"
        hist.y_label     = "Events"
        hist.format      = 'png'
        hist.saveAs      = self.output+"/hist_DNNscore_"+self.date
        hist.ATLASlabel       = 'top left'
        hist.ATLASlabelStatus = 'Simulation Internal'
        hist.numLegendColumns = 1
        hist.extra_text.Add(self.processlabel_args['text'],**self.processlabel_args)

        hist.initialize()

        # Train
        hist.Add(top_train_scores,name='score_top_train',linecolor='r',color='r',
                 draw='step',label='Large-R Jet (top) Train')
        hist.Add(tbar_train_scores,name='score_tbar_train',linecolor='b',color='b',
                 draw='step',label='Large-R Jet (anti-top) Train')

        # Test
        hist.Add(top_test_scores,name='score_top_test',linecolor='r',color='r',
                 draw='stepfilled',label='Large-R Jet (top) Test',alpha=0.5,linewidth=0)
        hist.Add(tbar_test_scores,name='score_tbar_test',linecolor='b',color='b',
                 draw='stepfilled',label='Large-R Jet (anti-top) Test',alpha=0.5,linewidth=0)

        p = hist.execute()
        hist.savefig()

        return


    def getROC(self):
        """Get the ROC curve"""
        self.fpr, self.tpr, _  = roc_curve(self.test_data['Y'], self.score)

        return


    def plot_ROC(self):
        """Plot the ROC curve"""
        self.getROC()

        fig,ax = plt.subplots()

        roc_auc = auc(self.fpr, self.tpr)
        ax.plot(self.fpr, self.tpr, label='ROC curve (AUC = %0.2f)' % roc_auc,lw=2)
        ax.plot([0, 1], [0, 1], 'k--',label='No Discrimination',lw=2)

        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.2])

        ax.set_xlabel('False Positive Rate',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(ax.get_xticks(),fontsize=22)
        ax.set_ylabel('True Positive Rate',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+list( ax.get_yticks()[1:-1] )+[''],fontsize=22)

        text_args = {'fontsize':18,'ha':'left','va':'top','transform':ax.transAxes}
        ax.text(0.03,0.97,r"\textbf{\textit{ATLAS}} Simulation Internal",**text_args)
        ax.text(0.03,0.90,r"$\sqrt{\text{s}}$ = 13 TeV",**text_args)

        leg = ax.legend(loc=4,numpoints=1,fontsize=18,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig(self.output+'/roc_curve_'+self.date+'.png',format='png',bbox_inches='tight',dpi=200)
        plt.close()

        return



## THE END ##