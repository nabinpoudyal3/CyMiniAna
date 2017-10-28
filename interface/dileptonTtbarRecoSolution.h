#ifndef DILEPTONTTBARRECOSOLUTION_H
#define DILEPTONTTBARRECOSOLUTION_H

#include <map>
#include "diHiggs/CyMiniAna/interface/tools.h"
#include "diHiggs/CyMiniAna/interface/physicsObjects.h"



class dileptonTtbarRecoSolution{

  public:

    /// Empty constructor
    dileptonTtbarRecoSolution();

    /// Destructor
    ~dileptonTtbarRecoSolution(){}


    /// Add a vector of solutions
    void Add(const std::vector<TtbarDilepton>& solutions);

    ///  Add a solution
    void Add(const TtbarDilepton& solution);


    /// Number of all solutions
    size_t numberOfSolutions() const {return v_solution_.size();}
    size_t numberOfSolutions(const int nBtags) const;

    /// Access from all solutions the one selected with solutionNumber, ranked by decreasing specific weight
    const TtbarDilepton solution(const TtbarDilepton::WeightType weightType, const size_t solutionNumber) const;

    /// Access from solutions with N b-tags the one selected with solutionNumber, ranked by decreasing specific weight
    TtbarDilepton solution(const TtbarDilepton::WeightType weightType,
                           const size_t solutionNumber, const int nBtags) const;

  private:

    /// Insert solution index for specific weight in vector for all solutions, ranked by weight
    void insertIndex(const size_t solutionIndex,
                     const double weight, std::vector<size_t>& v_index) const;

    /// Insert solution index for specific weight in vector for b-tag categorised solutions, ranked by weight
    void insertIndexByCategory(const std::vector<size_t>& v_solutionIndex,
                               const double weight, std::vector<size_t>& v_solutionIndexByCategory) const;



    /// Vector containing all solutions
    std::vector<TtbarDilepton> v_solution_;

    /// Vector containing indices of the solutions with 2 b-tags stored in v_solution_
    std::vector<size_t> v_solutionTwoBtags_;

    /// Vector containing indices of the solutions with 1 b-tag stored in v_solution_
    std::vector<size_t> v_solutionOneBtag_;

    /// Vector containing indices of the solutions with 0 b-tags stored in v_solution_
    std::vector<size_t> v_solutionNoBtags_;


    /// Map associating specific weight type to vector containing indices of all solutions, ordered for this weight
    std::map<TtbarDilepton::WeightType, std::vector<size_t>> m_weightIndex_;

    /// Map associating specific weight type to vector containing indices of 2 b-tag solutions, ordered for this weight
    std::map<TtbarDilepton::WeightType, std::vector<size_t>> m_weightIndexTwoBtags_;

    /// Map associating specific weight type to vector containing indices of 1 b-tag solutions, ordered for this weight
    std::map<TtbarDilepton::WeightType, std::vector<size_t>> m_weightIndexOneBtag_;

    /// Map associating specific weight type to vector containing indices of 0 b-tag solutions, ordered for this weight
    std::map<TtbarDilepton::WeightType, std::vector<size_t>> m_weightIndexNoBtags_;
};

#endif

