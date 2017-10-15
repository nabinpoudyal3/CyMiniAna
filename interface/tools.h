#ifndef TOOLS_H_
#define TOOLS_H_

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <assert.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TList.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

namespace cma{

    /* read a file and dump contents into vector of strings
       each line is a new element in the vector -- default comment = '#' */
    void read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment="#");
    std::ifstream open_file(const std::string &filename);
    void check_file(const std::ifstream& file, const std::string& fname);
    void check_file(const std::string & filename);

    /* Split a string with some delimiter (comma) */
    void split(const std::string &s, char delim, std::vector<std::string> &elems);

    /* Get the list of TTrees in a file */
    void getListOfKeys( TFile* file, std::vector<std::string> &fileKeys );

    /* calculate values for normalizing monte carlo samples */
    void getSampleWeights( std::string &sw_files, 
                           std::map<unsigned int,float> &m_mapXSection,
                           std::map<unsigned int,float> &m_mapKFactor,
                           std::map<unsigned int,float> &m_mapAMI );

    /* Convert string to boolean */
    bool str2bool( const std::string value );

    /* Convert vector of strings into a string of comma-separated elements */
    std::string vectorToStr( const std::vector<std::string> &vec );

    /* Calculate the median of a vector */
//    double median( std::vector<double> scores );

    template<typename T>
    T median( std::vector<T> scores ) {
        /* Calculate the median for a vector of values */
        T med;
        std::size_t size = scores.size();
        std::sort(scores.begin(), scores.end());
        if (size%2 == 0){ med = (scores[size / 2 - 1] + scores[size / 2]) / 2; }
        else{ med = scores[size / 2]; }
        return med;
    };

    /* DeltaR matching of TLorentzVectors (default deltaR=0.75) */
    bool deltaRMatch( TLorentzVector &particle1, TLorentzVector &particle2, double deltaR=0.75 );

    // debug message handling
    extern std::string m_debugLevel;
    void setVerboseLevel(const std::string& verboseLevel);
    void DEBUG(const std::string& message);
    void INFO(const std::string& message);
    void WARNING(const std::string& message);
    void ERROR(const std::string& message);
    void HELP(const std::string& runExecutable="run");
    std::map<std::string,unsigned int> verboseMap();
    void verbose(const std::string level, const std::string& message);
}

#endif
