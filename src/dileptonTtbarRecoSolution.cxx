/*
Created:        --
Last Updated:   22 October  2017

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University

-----

 Methods for dileptonTtbarReco
  Class for dilepton ttbar reconstruction.
  Imported from
   https://gitlab.cern.ch/cms-desy-top/TopAnalysis/blob/master/
           Configuration/analysis/common/src/KinematicReconstructionSolution.cc
  on 15 October 2017

*/
#include "cms-ttbarAC/CyMiniAna/interface/dileptonTtbarRecoSolution.h"


dileptonTtbarRecoSolution::dileptonTtbarRecoSolution() {}


void dileptonTtbarRecoSolution::Add(const std::vector<ttbarDilepton>& solutions) {
    /* Add solution to vector of solutions */
    for(const auto& solution : solutions) 
        this->Add(solution);

    return;
}


void dileptonTtbarRecoSolution::Add(const ttbarDilepton& solution) {
    // Initialise all maps of weight-ordered solution indices, in first solution only
    if(!m_weightIndex_.size()){
        for(const auto weightTypeWeight : solution.mapOfWeights){
            const auto weightType = weightTypeWeight.first;
            m_weightIndex_[weightType] = std::vector<size_t>();
            m_weightIndexTwoBtags_[weightType] = std::vector<size_t>();
            m_weightIndexOneBtag_[weightType]  = std::vector<size_t>();
            m_weightIndexNoBtags_[weightType]  = std::vector<size_t>();
        }
    }

    // Set pointers to the specific b-tag multiplicity category
    std::vector<size_t>* v_solutionByCategory(0);
    std::map<ttbarDilepton::WeightType, std::vector<size_t>>* m_weightIndexByCategory(0);
    const int numberOfBtags = solution.numberOfBtags();
    if(numberOfBtags == 2) {
        v_solutionByCategory    = &v_solutionTwoBtags_;
        m_weightIndexByCategory = &m_weightIndexTwoBtags_;
    }
    else if(numberOfBtags == 1) {
        v_solutionByCategory    = &v_solutionOneBtag_;
        m_weightIndexByCategory = &m_weightIndexOneBtag_;
    }
    else if(numberOfBtags == 0) {
        v_solutionByCategory    = &v_solutionNoBtags_;
        m_weightIndexByCategory = &m_weightIndexNoBtags_;
    }
    else{
        cma::ERROR("SOLUTION : Invalid number of b-tags: "+std::to_string(numberOfBtags));
        exit(731);
    }

    // Add solution to all solutions, and to specific b-tag multiplicity category
    v_solution_.push_back(solution);
    const size_t solutionIndex = std::distance(v_solution_.begin(), v_solution_.end()) - 1;
    v_solutionByCategory->push_back(solutionIndex);

    // Fill for each weight type the indices of solutions, ordered for the weight
    for(const auto weightTypeWeight : solution.weightMap()){
        const ttbarDilepton::WeightType weightType = weightTypeWeight.first;
        const double weight = weightTypeWeight.second;
        insertIndex(solutionIndex, weight, m_weightIndex_.at(weightType));
        insertIndexByCategory(*v_solutionByCategory, weight, m_weightIndexByCategory->at(weightType));
    }

    return;
}


void dileptonTtbarRecoSolution::insertIndex(const size_t solutionIndex,
                                             const double weight, 
                                             std::vector<size_t>& v_index) const {
    if(!v_index.size()){
        v_index.push_back(solutionIndex);
        return;
    }

    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_index.begin(); i_index != v_index.end(); ++i_index){
        const double& iWeight = v_solution_.at(*i_index).weight();
        if(iWeight < weight){
            v_index.insert(i_index, solutionIndex);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_index.push_back(solutionIndex);

    return;
}


void dileptonTtbarRecoSolution::insertIndexByCategory(const std::vector<size_t>& v_solutionIndex,
                                                       const double weight, 
                                                       std::vector<size_t>& v_solutionIndexByCategory) const {
    const size_t solutionIndexByCategory = std::distance(v_solutionIndex.begin(), v_solutionIndex.end()) - 1;

    if(!v_solutionIndexByCategory.size()){
        v_solutionIndexByCategory.push_back(solutionIndexByCategory);
        return;
    }

    bool isInserted(false);
    for(std::vector<size_t>::iterator i_index = v_solutionIndexByCategory.begin(); i_index != v_solutionIndexByCategory.end(); ++i_index){
        const double& iWeight = v_solution_.at(v_solutionIndex.at(*i_index)).weight();
        if(iWeight < weight){
            v_solutionIndexByCategory.insert(i_index, solutionIndexByCategory);
            isInserted = true;
            break;
        }
    }
    if(!isInserted) v_solutionIndexByCategory.push_back(solutionIndexByCategory);

    return;
}


size_t dileptonTtbarRecoSolution::numberOfSolutions(const int nBtags) const {
    /* Return the number of solutions for particular number of b-tags */
    size_t n_btag_soln = 0;
    if (nBtags==0){
        n_btag_soln = v_solutionNoBtags_.size();
    }
    else if (nBtags==1){
        n_btag_soln = v_solutionOneBtag_.size();
    }
    else if (nBtags==2){
        n_btag_soln = v_solutionTwoBtags_.size();
    }
    else {
        cma::ERROR("SOLUTION : Invalid number of b-tags: "+std::to_string(nBtags));
        exit(731);
    }

    return n_btag_soln;
}


const ttbarDilepton& dileptonTtbarRecoSolution::solution(const ttbarDilepton::WeightType weightType,
                                                          const size_t solutionNumber) const {
     /* Return particular solution */
    const size_t index = m_weightIndex_.at(weightType).at(solutionNumber);
    return v_solution_.at(index);
}


const ttbarDilepton& dileptonTtbarRecoSolution::solution(const ttbarDilepton::WeightType weightType,
                                                          const size_t solutionNumber,
                                                          const int nBtags) const {
     /* Return particular solution for specific number of b-tags */
    ttbarDilepton n_btag_soln;
    if (nBtags==0){
        const size_t index = m_weightIndexNoBtags_.at(weightType).at(solutionNumber);
        n_btag_soln = v_solution_.at(v_solutionNoBtags_.at(index));
    }
    else if (nBtags==1){
        const size_t index = m_weightIndexOneBtags_.at(weightType).at(solutionNumber);
        n_btag_soln = v_solution_.at(v_solutionOneBtags_.at(index));
    }
    else if (nBtags==2){
        const size_t index = m_weightIndexTwoBtags_.at(weightType).at(solutionNumber);
        n_btag_soln = v_solution_.at(v_solutionTwoBtags_.at(index));
    }
    else {
        cma::ERROR("SOLUTION : Invalid number of b-tags: "+std::to_string(nBtags));
        exit(731);
    }

    return n_btag_soln;
}

// THE END