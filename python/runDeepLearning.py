"""
Created:        12 November  2016
Last Updated:   11 December  2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109
-----

Script for running the deep learning implementation
Top vs Anti-top identification
"""
import os
import sys
from deepLearning import DeepLearning
from argparse import ArgumentParser
from time import strftime,localtime

date   = strftime("%d%b", localtime())
cmaDir = os.path.expandvars('$CYMINIANADIR')
numbers = {1:'one',2:'two',3:'three',4:'four',5:'five',6:'six',60:'sixty'}

## Argument parser to set configuration options ##
parser = ArgumentParser(description="Deep Learning")

parser.add_argument('--buildNN', action='store_true',
                    dest='buildNN',default=False,
                    help='Boolean: build DNN')
parser.add_argument('--loadNN', action='store_true',
                    dest='loadNN',default=False,
                    help='Boolean: load DNN')
parser.add_argument('--hep_data', action='store',
                    dest='hep_data',default=cmaDir+'/data/ttbar_allhad_ttbarXS_77.json',
                    help='Physics data to load for NN')
parser.add_argument('--dnn_data', action='store',
                    dest='dnn_data',default=cmaDir+'/data/ttbar_allhad_ttbarXS_77_'+date,
                    help='Base of filename for loading .json and .h5 NN model')
parser.add_argument('--nHiddenLayers',action='store',type=int,
                    dest='nHiddenLayers',default=2,
                    help='Number of hidden layers in Neural Network')
parser.add_argument('--verbose', action='store_true',
                    dest='verbose',default=False,
                    help='Boolean: Print model output')
parser.add_argument('--output', action='store',
                    dest='output',default=cmaDir+'/data/DNN/one_layer/zprime_training_'+date+'/ttbar_had_testing/',
                    help='Output directory')
config = parser.parse_args()




print " ------------------------------ "
print " *  Deep Learning with Keras  * "
print " ------------------------------ "
print " Discrimination between hadronic "
print " decays of top/antitop quarks "
print


dnn = DeepLearning()
dnn.hep_data  = config.hep_data
dnn.dnn_data  = config.dnn_data
nHiddenLayers = numbers[config.nHiddenLayers]
dnn.output    = 'data/DNN/{0}_layers/{1}'.format(nHiddenLayers,config.hep_data.split('/')[-1].split('.')[0])  # config.output

if not os.path.isdir(dnn.output):
    print " WARNING : '{0}' does not exist ".format(dnn.output)
    print "           Creating the directory. "
    os.system( 'mkdir -p {0}'.format(dnn.output) )

# load hep data (physics data -- .json file). Always need this for testing/training
dnn.getHEPData()
dnn.plot_features()  # sanity checks
dnn.save_features()  # save the features for loading with lwtnn later

if config.buildNN:
    outputFile   = open(dnn.output+'/ABOUT.txt','w')
    outputFile.write(" * NN Setup * \n")

    NN_parameters  = ['nb_epoch','batch_size','loss','optimizer','metrics','nHiddenLayers','output_dim','input_dim']
    print " > Build the NN"
    # set properties of the NN
    dnn.nb_epoch   = 10
    dnn.batch_size = 32
    dnn.loss       = 'binary_crossentropy'
    dnn.optimizer  = 'adam'
    dnn.metrics    = ['accuracy']
    dnn.nHiddenLayers = config.nHiddenLayers

    dnn.buildNN()

    outputFile.write(" NN Summary: \n")
    dnn.model.summary()
    outputFile.write("\n NN parameters: \n")
    for NN_parameter in NN_parameters:
        outputFile.write( NN_parameter+": "+str(getattr(dnn,NN_parameter))+"\n" )

    outputFile.write( "\n NN Features: \n" )
    for feature in dnn.features:
        outputFile.write("  >> "+feature+"\n" )

if config.loadNN:
    print " > Load NN model from disk"
    dnn.loadModel()

if not config.loadNN and not config.buildNN:
    print
    print " No configuration set "
    print " Please use the arguments 'loadNN' or 'buildNN' to define workflow "
    results.help()


print " > Get NN score "
score = dnn.getScore()
dnn.plot_score()

print " > Obtain & plot the ROC curve"
dnn.plot_ROC()

if config.buildNN:
    print " > Save model. "
    dnn.saveModel()

