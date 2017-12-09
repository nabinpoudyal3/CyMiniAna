"""
Created:        11 November  2016
Last Updated:   11 December  2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Class for performing deep learning to classify 
large-R jets as either a top or anti-top quark.

https://keras.io/
"""
import json
import util
import datetime
import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
import hepPlotterLabels as hpl
import hepPlotterTools as hpt
from hepPlotter import HepPlotter
import keras
from keras.models import Sequential,model_from_json
from keras.layers import Dense, Activation
from keras.utils.visualize_util import plot
from keras.callbacks import EarlyStopping
from sklearn.model_selection import train_test_split,StratifiedKFold
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
        hpl_text_dicts  = hpl.text_dicts()
        self.text_dicts = hpl_text_dicts['dnn']
        self.processlabel_args = hpl_text_dicts['samples']

        self.text_args  = {'fontsize':18,'ha':'left','va':'top','transform':None}
        self.verbose    = True
        self.history    = None
        self.date       = datetime.date.today().strftime('%d%b%Y')
        self.hep_data   = None
        self.dnn_data   = None
        self.df         = None # set later
        self.model      = None # set later
        self.fpr        = None # set later
        self.tpr        = None # set later
        self.accuracy   = {'mean':0,'std':0}   # set later
        self.metadata   = {}   # set later
        self.train_data = {}   # set later
        self.test_data  = {}   # set later
        self.rejections = []   # set later
        self.rejection  = {}   # set later
        self.percentile = 75
        self.output     = 'data/DNN/' # set later
        self.dnn_name   = "ttbar_DNN"  # name to access in lwtnn ('variables.json')
        self.image_format = 'png'
        self.test_scores  = []   # set later
        self.train_scores = []   # set later
        self.verbose_level = "INFO"

        ## NN parameters -- set by config file
        self.features   = []
        self.batch_size = 32
        self.nb_epoch   = 1        
        self.loss       = 'binary_crossentropy' # preferred for binary classification
        self.optimizer  = 'adam' # old: 'sgd'; change to 'adam' per Riccardo di Sipio
        self.metrics    = ['accuracy']
        self.init       = 'normal' #'uniform'
        self.input_dim  = len(self.features)
        self.nHiddenLayers = 1
        self.nNodes        = [5 for _ in range(self.nHiddenLayers)]
        self.activation    = 'elu'
        self.elu_alpha     = 1.0
        self.kfold_splits  = 2
        self.earlystopping = {'monitor':'loss','min_delta':0.0001,
                              'patience':5,'verbose':self.verbose,'mode':'auto'}

        return



    def initialize(self):
        """Initialize a few parameters that are required before running, but after __init__()"""
        self.msg_svc       = util.VERBOSE()
        self.msg_svc.level = self.verbose_level

        return

    def getHEPData(self):
        """
        Load the data for NN.  
        This uses data created by root2keras.py where there is metadata 
        saved to the json output (to reproduce results later &
        understand how they were produced) <- remove these extra keys.
        """
        self.features2plot = self.features+['pt','eta','reco_m_ttbar','truth_m_ttbar']

        data  = json.load( open(self.hep_data,'r') )
        fdata = dict(  (key,value) for key,value in data.items() if key in self.features2plot+['target'] )
        self.metadata = data['metadata']

        # -- convert to DataFrame for easier slicing 
        self.df = pd.DataFrame( fdata )

        ## custom modification for tjet_pt -> normalize by the large-R jet pT; max value for n_tjets = 3
        if 'tjet_0_pt' in self.features:
            self.df['tjet_0_pt'] = self.df['tjet_0_pt']/self.df['pt']
            if 'tjet_1_pt' in self.features:
                self.df['tjet_1_pt'] = self.df['tjet_1_pt']/self.df['pt']
            if 'tjet_2_pt' in self.features:
                self.df['tjet_2_pt'] = self.df['tjet_2_pt']/self.df['pt']

        self.test_data  = {'X':[],'Y':[]}
        self.train_data = {'X':[],'Y':[]}
        # -- split the dataset into train and test portions -- obsolete
        #    now handled by k-fold cross validation below.
        #    This is here for reference, if needed
        # X_data  = self.df[self.features].values
        # Y_data  = self.df['target'].values
        # X_train, X_test, y_train, y_test = train_test_split(X_data,Y_data,test_size=0.4)

        return



    def buildNN(self):
        """Initialize the NN"""

        if not self.train_data and not self.test_data:
            self.getHEPData()
        if len(self.nNodes)==1 and self.nHiddenLayers>1:
            nodes_per_layer = self.nNodes[0]
            self.nNodes = [nodes_per_layer for _ in range(self.nHiddenLayers)]

        # change the activation function (from Riccardo di Sipio)
        elu = keras.layers.advanced_activations.ELU( alpha=self.elu_alpha )

        activations = {'elu':elu,'relu':Activation('relu')}

        # define early stopping callback
        earlystop = EarlyStopping(**self.earlystopping)
        callbacks_list = [earlystop]

        # Declare the model
        self.model = Sequential() # The Keras Sequential model is a linear stack of layers.

        # Add 1st layer
        self.model.add( Dense( int(self.nNodes[0]), input_dim=self.input_dim, init=self.init) )
        self.model.add( activations[self.activation] )  # old: Activation("relu") )

        # Add middle layer(s) [identical setup for now]
        # Number of middle layers = hidden layers - 1 (initial layer above)
        for h in range(self.nHiddenLayers-1):
            self.model.add( Dense( int(self.nNodes[h+1]), init=self.init) )
            self.model.add( activations[self.activation] )   # old: Activation("relu") )

        # Add the output
        self.model.add( Dense(1,init=self.init) )
        self.model.add( Activation("sigmoid") )

        # Build the model
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)


        ## Setup for training the model using k-fold cross-validation ##

        # define k-fold cross validation test
        self.test_data  = {'X':[],'Y':[]}
        self.train_data = {'X':[],'Y':[]}
        self.test_scores = []
        self.train_scores = []
        cvscores   = []

        X = self.df[self.features].values
        Y = self.df['target'].values

        kfold = StratifiedKFold(n_splits=self.kfold_splits, shuffle=True, random_state=seed)
        kfold.get_n_splits(X,Y)

        for ind,(train,test) in enumerate(kfold.split(X,Y)):
            self.msg_svc.INFO("DL :   - Fitting K-Fold {0}".format(ind))

            self.test_data['X'].append(X[test])
            self.test_data['Y'].append(Y[test])
            self.train_data['X'].append(X[train])
            self.train_data['Y'].append(Y[train])

            md_info = self.model.fit(X[train],Y[train],nb_epoch=self.nb_epoch,
                                     callbacks=callbacks_list,batch_size=self.batch_size,verbose=self.verbose)

            # evaluate the model
            self.msg_svc.INFO("DL :     + Evaluate the model: ")
            scores = self.model.evaluate(X[test], Y[test],verbose=self.verbose,batch_size=self.batch_size)
            cvscores.append(scores[1] * 100)
            self.msg_SVc.INFO("DL :       {0}: {1:.2f}%".format(self.model.metrics_names[1], scores[1]*100))

            test_scores  = deepcopy( self.getScore(X[test]) )
            self.test_scores.append( test_scores )
            train_scores = deepcopy( self.getScore(X[train]) )
            self.train_scores.append( train_scores )

            self.plot_model_loss(md_info,kfold=ind,val_loss=scores[0])

        self.msg_svc.INFO("DL :   Finished K-Fold cross-validation: ")
        self.accuracy = {'mean':np.mean(cvscores),'std':np.std(cvscores)}
        self.msg_svc.INFO("DL :   - Accuracy: {0:.2f}% (+/- {1:.2f}%)".format(np.mean(cvscores), np.std(cvscores)))

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
        plot(self.model,to_file=output+'_model.eps',show_shapes=True)

        return



    def loadModel(self):
        """Load existing model to make plots or predictions"""
        json_file  = open(self.dnn_data+"_model.json", 'r')
        model_json = json_file.read()

        self.model = model_from_json(model_json)
        self.model.load_weights(self.dnn_data+"_weights.h5")

        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)

        return



    def getScore(self,data):
        """Return the score from a test sample"""
        self.msg_svc.DEBUG("DL : Get the DNN score")
        score = self.model.predict( data, batch_size=self.batch_size )

        return score


    def plot_features(self):
        """Plot the features"""
        self.msg_svc.INFO("DL : Plotting features comparing top quarks and anti-quarks. ")
        top  = self.df.loc[self.df['target'] == self.metadata['t_target']]
        tbar = self.df.loc[self.df['target'] == self.metadata['tbar_target']]

        filename = self.metadata['file'].split('/')[-1].split('.')[0].rstrip('\n')

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
            hist.format      = self.image_format
            hist.saveAs      = self.output+"/hist_"+feature+"_"+self.date
            hist.ATLASlabel       = 'top left'
            hist.ATLASlabelStatus = 'Simulation Internal'
            hist.numLegendColumns = 1
            hist.extra_text.Add(self.processlabel_args[filename]['label'],coords=[0.03,0.80])

            hist.initialize()

            multiply=1.
            if feature.endswith('_m_ttbar') or feature=='pt':
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
        namelabels = [r"t correlations",r"$\bar{\text{t}}$ correlations"]
        fontProperties = {'family':'sans-serif'}
        opts = {'cmap': plt.get_cmap("bwr"), 'vmin': -1, 'vmax': +1}

        for c,corrmat in enumerate([corrmat_df_top,corrmat_df_tbar]):

            fig,ax = plt.subplots()

            # hide the upper part of the triangle
            #mask = np.zeros_like(corrmat, dtype=np.bool)    # return array of zeros with same shape as corrmat
            #mask[np.tril_indices_from(mask)] = True
            #corrmat_mask  = np.ma.array(corrmat, mask=mask) 

            heatmap1 = ax.pcolor(corrmat, **opts)
            cbar     = plt.colorbar(heatmap1, ax=ax)

            cbar.ax.set_yticklabels( [i.get_text().strip('$') for i in cbar.ax.get_yticklabels()], **fontProperties )

            labels = corrmat.columns.values
            labels = [i.replace('_','\_') for i in labels]
            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_xticklabels(labels, fontProperties, fontsize=18, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(labels, fontProperties, fontsize=18, minor=False)

            text_args = {'fontsize':16,'ha':'left','va':'bottom','transform':ax.transAxes}

            ## ATLAS Label + Signal name
            ax.text(0.02,1.00,r"\textbf{\textit{ATLAS}} Simulation Internal",**text_args)
            ax.text(0.03,0.93,"{0}, {1}".format(self.processlabel_args[filename]['label'],namelabels[c]),**text_args)

            ## Energy Label
            text_args['ha'] = 'right'
            ax.text(0.99,1.00,r"$\sqrt{\text{s}}$ = 13 TeV",**text_args)

            plt.savefig(self.output+"/correlations_{0}_{1}.{2}".format(names[c],self.date,self.image_format),
                        format=self.image_format,dpi=300,bbox_inches='tight')
            plt.close()

        return


    def save_features(self):
        """
        Save the features to a json file to load in the lwtnn later
        Hard-coded scale & offset; must change later if necessary
        """
        text = """  {
    "inputs": ["""

        for fe,feature in enumerate(self.features):
            comma = "," if fe!=self.features else ""
            tmp = """
      {"name": "%(feature)s",
       "scale":  1,
       "offset": 0}%(comma)s""" % {'feature':feature,'comma':comma}
            text += tmp
        text += "],"
        text += """
    "class_labels": ["%(name)s"],
    "keras_version": "%(version)s",
    "miscellaneous": {"elu_alpha":%(alpha)1.1f}
  }
""" % {'version':keras.__version__,'name':self.dnn_name,'alpha':self.elu_alpha}

        varsFileName = self.output+'/variables.json'
        varsFile     = open(varsFileName,'w')
        varsFile.write(text)

        return


    def plot_score(self):
        """Plot the features"""

        betterColors = hpt.betterColors()['linecolors']
        filename     = self.metadata['file'].split('/')[-1].split('.')[0].rstrip('\n')

        # Plot all k-fold cross-validation results
        for i,(train_X,train_Y,test_X,test_Y) in enumerate(zip(self.train_data['X'],self.train_data['Y'],self.test_data['X'],self.test_data['Y'])):

            hist = HepPlotter("histogram",1)

            hist.ratio_plot  = True
            hist.y_ratio_label = "Test/Train"
            hist.normed      = True
            hist.binning     = [0.05*j for j in range(21)]
            hist.stacked     = False
            hist.logplot     = False
            hist.x_label     = "DNN Score"
            hist.y_label     = "Events"
            hist.format      = self.image_format
            hist.label_size  = 14
            hist.saveAs      = self.output+"/hist_DNNscore_kfold{0}_{1}".format(i,self.date)
            hist.ATLASlabel       = 'top left'
            hist.ATLASlabelStatus = 'Simulation Internal'
            hist.numLegendColumns = 1
            hist.extra_text.Add(self.processlabel_args[filename]['label'],coords=[0.03,0.80],fontsize=14)

            hist.initialize()

            top_train_scores  = self.train_scores[i][ train_Y==1 ]
            tbar_train_scores = self.train_scores[i][ train_Y==0 ]

            top_test_scores   = self.test_scores[i][ test_Y==1 ]
            tbar_test_scores  = self.test_scores[i][ test_Y==0 ]

            ## Train
            index = i*2
            top_color  = 'r' #betterColors[index]
            tbar_color = 'b' #betterColors[index+1]
            hist.Add(top_train_scores,name='score_top_train_'+str(i),linecolor=top_color,color=top_color,linewidth=2,
                     draw='step',label='Large-R Jet (top) Train '+str(i),ratio_den=True,ratio_num=False,ratio_partner='score_top_test_'+str(i))
            hist.Add(tbar_train_scores,name='score_tbar_train_'+str(i),linecolor=tbar_color,color=tbar_color,linewidth=2,
                     draw='step',label='Large-R Jet (anti-top) Train '+str(i),ratio_den=True,ratio_num=False,ratio_partner='score_tbar_test_'+str(i))

            ## Test
            hist.Add(top_test_scores,name='score_top_test_'+str(i),linecolor=top_color,color=top_color,
                     draw='stepfilled',label='Large-R Jet (top) Test '+str(i),alpha=0.5,linewidth=0,ratio_num=True,ratio_den=False,ratio_partner='score_top_train_'+str(i))
            hist.Add(tbar_test_scores,name='score_tbar_test_'+str(i),linecolor=tbar_color,color=tbar_color,
                     draw='stepfilled',label='Large-R Jet (anti-top) Test '+str(i),alpha=0.5,linewidth=0,ratio_num=True,ratio_den=False,ratio_partner='score_tbar_train_'+str(i))

            p = hist.execute()
            hist.savefig()

            ## Calculation of the rejection
            ## use percentile (set above) to calculate at specific efficiency
            eff_value  = np.percentile(top_test_scores,self.percentile)
            tbar_wrong = tbar_test_scores[tbar_test_scores>=eff_value]
            rejection  = len(tbar_wrong)*1.0 / len(tbar_test_scores)
            self.rejections.append(rejection)

        self.rejection = {'mean':np.mean(self.rejections),'std':np.std(self.rejections)}

        return


    def getROC(self):
        """Get the ROC curve"""
        self.fpr = []  # false positive rate
        self.tpr = []  # true positive rate
        for s,score in enumerate(self.test_scores):
            fpr,tpr,_  = roc_curve(self.test_data['Y'][s],score)
            self.fpr.append(fpr)
            self.tpr.append(tpr)

        return


    def plot_ROC(self):
        """Plot the ROC curve"""
        self.getROC()

        fig,ax = plt.subplots()

        # Draw all of the ROC curves from the K-fold cross-validation
        ax.plot([0, 1], [0, 1], ls='--',label='No Discrimination',lw=2,c='gray')
        for ft,(fpr,tpr) in enumerate(zip(self.fpr,self.tpr)):
            roc_auc = auc(fpr,tpr)
            ax.plot(fpr, tpr, label='K-fold {0} (AUC = {1:.2f})'.format(ft,roc_auc),lw=2)

        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.5])

        ax.set_xlabel(r'$\epsilon$(anti-top)',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(ax.get_xticks(),fontsize=22)
        ax.set_ylabel(r'$\epsilon$(top)',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+list( ax.get_yticks()[1:-1] )+[''],fontsize=22)

        filename = self.metadata['file'].split('/')[-1].split('.')[0].rstrip('\n')
        self.text_args['transform'] = ax.transAxes
        ax.text(0.03,0.97,r"\textbf{\textit{ATLAS}} Simulation Internal",**self.text_args)
        ax.text(0.03,0.90,r"$\sqrt{\text{s}}$ = 13 TeV",**self.text_args)
        ax.text(0.03,0.82,self.processlabel_args[filename]['label'],**self.text_args)
        ax.text(0.03,0.75,r"Accuracy = {0:.2f}$\pm${1:.2f}".format(self.accuracy['mean'],self.accuracy['std']),**self.text_args)
        #ax.text(0.03,0.68,r"Rejection ({0}\% $\epsilon$) = {1:.2f}".format(self.percentile,self.rejection['mean']),**self.text_args)
        ax.axhline(y=1,lw=1,c='lightgray',ls='--')

        leg = ax.legend(loc=4,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig(self.output+'/roc_curve_{0}.{1}'.format(self.date,self.image_format),
                    format=self.image_format,bbox_inches='tight',dpi=300)
        plt.close()

        return


    def plot_model_loss(self,model_history,kfold=0,val_loss=0.0):
        """Plot loss as a function of epoch for model"""
        fig,ax = plt.subplots()

        loss = model_history.history['loss']
        ax.plot(range(1,len(loss)+1),loss,    label='Training',  color='r')

        ax.set_xlabel('Epoch',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(ax.get_xticks(),fontsize=22)
        ax.set_ylabel('Loss',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+list( ax.get_yticks()[1:-1] )+[''],fontsize=22)

        filename = self.metadata['file'].split('/')[-1].split('.')[0].rstrip('\n')
        self.text_args['transform'] = ax.transAxes
        ax.text(0.03,0.97,r"\textbf{\textit{ATLAS}} Simulation Internal",**self.text_args)
        ax.text(0.03,0.90,r"$\sqrt{\text{s}}$ = 13 TeV",**self.text_args)
        ax.text(0.03,0.82,"{0}".format(self.processlabel_args[filename]['label'],kfold),**self.text_args)
        ax.text(0.03,0.76,"Validation Loss = {0}; K-fold {1}".format(val_loss,kfold),**self.text_args)

        leg = ax.legend(loc=1,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig(self.output+'/loss_epochs_{0}_{1}.{2}'.format(kfold,self.date,self.image_format),
                    format=self.image_format,bbox_inches='tight',dpi=200)
        plt.close()


## THE END ##

