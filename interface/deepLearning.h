#ifndef DEEPLEARNING_H
#define DEEPLEARNING_H

#include <string>
#include <map>
#include <vector>

#include "lwtnn/lwtnn/interface/LightweightNeuralNetwork.hh"
#include "lwtnn/lwtnn/interface/parse_json.hh"

#include "Analysis/CyMiniAna/interface/tools.h"
#include "Analysis/CyMiniAna/interface/configuration.h"
#include "Analysis/CyMiniAna/interface/physicsObjects.h"


class deepLearning {
  public:
    deepLearning( configuration& cmaConfig );

    ~deepLearning();

    void training(std::vector<Ljet>& ljets);
    void training(Ljet& ljet);

    void inference(std::vector<Ljet>& ljets);
    void inference(Ljet& ljet);

    void loadFeatures(const Ljet& ljet);

    std::map<std::string,double> predictions() const {return m_predictions;}
    double prediction() const {return m_DNN;}
    double prediction(const std::string& key) const;

    std::map<std::string,double> features() const {return m_features;}

  protected:

    configuration *m_config;

    lwt::LightweightNeuralNetwork* m_lwnn;       // LWTNN tool

    std::map<std::string, double> m_features;    // values for inputs to the DNN
    std::map<std::string,double> m_predictions;  // map of DNN predictions
    std::string m_dnnKey;                        // default key for accessing map of values
    float m_DNN;                                 // DNN prediction for one key
};

#endif

