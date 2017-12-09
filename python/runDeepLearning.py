"""
Created:        12 November  2016
Last Updated:    2 April     2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Script for running the deep learning implementation
Top vs Anti-top identification

To run:
$ python python/runDeepLearning.py share/mlconfig.txt
-- the second argument is the text file with configurations for the NN/setup
"""
import os
import sys
import json
from info import VERBOSE
from config import Config
from collections import Counter
from time import strftime,localtime
from deepLearning import DeepLearning


print
print " ------------------------------ "
print " *  Deep Learning with Keras  * "
print " ------------------------------ "
print " Discrimination between hadronic "
print " decays of top/antitop quarks "
print


date   = strftime("%d%b", localtime())
cmaDir = os.path.expandvars('$CYMINIANADIR')
vb     = VERBOSE()

## Set configuration options ##
config = Config(sys.argv[1])
vb.level = config.verbose_level


if not config.loadNN and not config.buildNN:
    vb.ERROR("RUN :  No configuration set ")
    vb.ERROR("RUN :  Please set the arguments 'loadNN' or 'buildNN' to define workflow ")
    vb.ERROR("RUN :  Exiting.")
    sys.exit(1)


## Setup features
featureKeys = json.load(open('share/features.json','r'))

featureKey = -1
for key in featureKeys.keys():
    if Counter(featureKeys[key])==Counter(config.features):
        featureKey = int(key)
        break
if featureKey<0:
    featureKey = max([int(i) for i in featureKeys.keys()])+1
    featureKeys[str(featureKey)] = config.features
    vb.INFO("RUN :  New features for NN ")
    with open('share/features.json','w') as outfile:
        json.dump(featureKeys,outfile)


## Setup Deep Learning class
dnn = DeepLearning()
dnn.hep_data  = config.hep_data
dnn.dnn_data  = config.dnn_data
dnn.verbose_level = config.verbose_level



## Set output directory
output_dir  = "nHiddenLayers{0}_".format(config.nHiddenLayers)
output_dir += "nNodes{0}_".format('-'.join(config.nNodes))
output_dir += "epoch{0}_".format(config.nb_epoch)
output_dir += "batch{0}_".format(config.batch_size)
output_dir += "kfold{0}_".format(config.kfold_splits)
output_dir += "activation-{0}_".format(config.activation)
output_dir += "featureKey{0}".format(featureKey)
hep_data_name = config.hep_data.split('/')[-1].split('.')[0]
if config.buildNN:
    dnn.output = config.output_path+"/{0}/{1}".format(output_dir,hep_data_name)
else:
    dnn.output = config.output_path+"/{0}/{1}/loadDNN".format(output_dir,hep_data_name)

if not os.path.isdir(dnn.output):
    vb.WARNING("RUN : '{0}' does not exist ".format(dnn.output))
    vb.WARNING("RUN :       Creating the directory. ")
    os.system( 'mkdir -p {0}'.format(dnn.output) )
else:
    vb.INFO("RUN :  Saving output to {0}".format(dnn.output))



## load hep data (physics data -- .json file). Always need this for testing/training
dnn.features = config.features
dnn.initialize()     # Setup some things as needed
dnn.getHEPData()     # load JSON file with HEP data
dnn.plot_features()  # sanity checks
dnn.save_features()  # save the features for loading with lwtnn later

if config.buildNN:

    vb.INFO("RUN :  > Build the NN")
    # set properties of the NN
    dnn.nb_epoch   = config.nb_epoch
    dnn.batch_size = config.batch_size
    dnn.loss       = config.loss
    dnn.optimizer  = config.optimizer
    dnn.metrics    = config.metrics
    dnn.init       = config.init
    dnn.input_dim  = len(config.features)
    dnn.nNodes     = config.nNodes
    dnn.nHiddenLayers = config.nHiddenLayers
    dnn.kfold_splits  = config.kfold_splits
    dnn.percentile    = config.percentile
    dnn.activation    = config.activation

    dnn.buildNN()
    dnn.model.summary()


    # Save information on the NN to a text file to reference later
    outputFile = open(dnn.output+'/ABOUT.txt','w')
    outputFile.write(" * NN Setup * \n")
    outputFile.write(" NN Summary: \n")
    outputFile.write("\n NN parameters: \n")

    NN_parameters = ['nb_epoch','batch_size','loss','optimizer','metrics','activation',
                     'nHiddenLayers','nNodes','input_dim','kfold_splits','percentile']
    for NN_parameter in NN_parameters:
        outputFile.write( NN_parameter+": "+str(getattr(dnn,NN_parameter))+"\n" )

    outputFile.write( "\n NN Features: \n" )
    for feature in dnn.features:
        outputFile.write("  >> "+feature+"\n" )

    outputFile.close()


if config.loadNN:
    vb.INFO("RUN :  > Load NN model from disk")
    dnn.loadModel()


vb.INFO("RUN :  > Get NN score ")
dnn.plot_score()

vb.INFO("RUN :  > Obtain & plot the ROC curve")
dnn.plot_ROC()

if config.buildNN:
    vb.INFO("RUN :  > Save model. ")
    dnn.saveModel()


## END ##
