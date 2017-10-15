/*
Created:        14 May       2016
Last Updated:   20 August    2016

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
University of Michigan, Ann Arbor, MI 48109

-----

Common tools needed

*/
#include "diHiggs/CyMiniAna/interface/tools.h"

namespace cma{


void check_file(const std::string & filename) {
    /* Check that file exists */
    std::ifstream f = open_file(filename);

    return;
}


void check_file(const std::ifstream& file, const std::string& fname){
    /* Check that file exists */

    if (!file) {
        cma::ERROR("TOOLS : File does not exist:       "+fname);
        cma::ERROR("TOOLS : Exiting. ");
        assert(file);
    }

    return;
}

std::ifstream open_file(const std::string &filename) {
    /* Open file */
    std::ifstream ifile(filename.c_str());
    check_file(ifile,filename);

    return ifile;
}


void read_file( const std::string &file_name, std::vector<std::string> &values, const std::string &comment ) {
    /* Read in a generic file and put it into a vector of strings */
    std::ifstream tmp_name = open_file(file_name);

    // open the file and put the data into a vector
    std::string line("");
    if (tmp_name.is_open()){

        while (std::getline(tmp_name, line)) {
            std::string newstring(line);

            // allow for comments
            std::size_t lineComment = line.find(comment);
            if (lineComment != std::string::npos) newstring = line.substr(0, lineComment);

            // remove all white spaces at the end of the string
            std::size_t space_pos = newstring.rfind(" ");
            while ( space_pos != std::string::npos && space_pos == newstring.size()-1 ) {
                newstring = newstring.substr(0, newstring.rfind(" "));
                space_pos = newstring.rfind(" ");
            }

            // ignore empty lines
            if(newstring.length()==0) continue;

            values.push_back(newstring); // put values into vector
        }

        tmp_name.close();
    }

    return;
}



void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }

    return;
}


void getListOfKeys( TFile* file, std::vector<std::string> &fileKeys ){
    /* Find the list of TTrees in a file */
    fileKeys.clear();

    TList *list = file->GetListOfKeys();
    TIter iter(list->MakeIterator());
    while(TObject* obj = iter()){
        TKey* key = (TKey*)obj;
        std::string keyname( key->GetName() );
        fileKeys.push_back(keyname);
    }

    return;
}


void getSampleWeights( std::string &sw_files, 
                       std::map<unsigned int,float> &m_mapXSection,
                       std::map<unsigned int,float> &m_mapKFactor,
                       std::map<unsigned int,float> &m_mapAMI ){
    /* Calculate XSection, KFactor, and sum of weights (AMI) */
    cma::INFO("TOOLS : Get sample weights (including sum of weights)");

    // get the absolute path (in case of batch job)
    char* cma_path = getenv("CYMINIANADIR");
    std::string cma_absPath("");
    if (cma_path==NULL){
        cma::WARNING("TOOLS : environment variable " );
        cma::WARNING("TOOLS :    'CYMINIANADIR' " );
        cma::WARNING("TOOLS : is not set.  Relative paths will be used " );
        cma::WARNING("TOOLS : This may cause problems submitting batch jobs." );
        cma_absPath = ".";
        cma::WARNING("TOOLS : path set to: "+cma_absPath );
    }
    else{
        cma_absPath = cma_path;
    }

    // Read in XSection data
    std::ifstream in = open_file("config/XSection-MC15-13TeV.data");

    for( ; !in.eof() ; ) {
        std::string line;
        if (!std::getline(in,line)) break;
        if (!line.empty() && line[0]!='#') {
            int dsid(-1);
            float xSect,kFact;
            std::string s_shower;

            std::istringstream istr(line);
            istr >> dsid >> xSect >> kFact >> s_shower;
            // ignoring everything after the KFactor

            m_mapXSection[dsid]  = xSect;
            m_mapKFactor[dsid]   = kFact;
        }
    }
    in.close();
    m_mapXSection[0]  = 1.0; // protection for Data
    m_mapKFactor[0]   = 1.0;

    // save sum of weights before derivations
    std::vector<std::string> sw_filenames;
    read_file( sw_files, sw_filenames, "#");

    float sum_of_weights(0.0);
    for (const auto& filename : sw_filenames) {

        auto file = TFile::Open(filename.c_str());
        if (!file || file->IsZombie()){
            cma::WARNING("TOOLS : File "+filename+" does not exist. Continuing.");
            continue;
        }

        TTreeReader sumWeights("sumWeights", file);
        TTreeReaderValue<float> totalEventsWeighted(sumWeights, "totalEventsWeighted");
        TTreeReaderValue<unsigned long long> totalEvents(sumWeights, "totalEvents");
        TTreeReaderValue<int>   dsid(sumWeights, "dsid");

        sum_of_weights = 0.0;
        while (sumWeights.Next()){
            if (*dsid >= 361020 && *dsid <= 361032)
                sum_of_weights = 1.0*(*totalEvents);    // different weighting for JZ*W samples
            else
                sum_of_weights = (*totalEventsWeighted);

            if (m_mapAMI.find((*dsid)) == m_mapAMI.end()){
                m_mapAMI[(*dsid)]  = sum_of_weights;
            }
            else{
                m_mapAMI[(*dsid)] += sum_of_weights;
            }
        }

        delete file;
        file = ((TFile *)0);
    } // end loop over files
    m_mapAMI[0] = 1.0; // protection for Data

    return;
}


bool str2bool( const std::string value ){
    /* Turn string into boolean */
    bool valueBoolean(false);

    if (value.compare("True")==0 || value.compare("true")==0 || value.compare("1")==0){
        valueBoolean = true;
    }
    else{
        valueBoolean = false;
    }

    return valueBoolean;
}


std::string vectorToStr( const std::vector<std::string> &vec ){
    std::string str_list;
    for( const std::string &str : vec)
        str_list += str + std::string(",");
    str_list.pop_back(); // remove last comma
    return str_list;
}


//template typename median<int>;
//template typename median<float>;
//template typename median<double>;
/*
double median(std::vector<double> scores) {
    // Calculate the median for a vector of values //
    double med;
    std::size_t size = scores.size();
    std::sort(scores.begin(), scores.end());

    if (size%2 == 0)
        med = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    else 
        med = scores[size / 2];

    return med;
}
*/

bool deltaRMatch( TLorentzVector &particle1, TLorentzVector &particle2, double deltaR ){
    /* Do the deltaR calculation (in one place) */
    bool isMatched(false);

    if (particle1.DeltaR( particle2 ) < deltaR){
        isMatched=true;
    }

    return isMatched;
}


std::string m_debugLevel = "SetMe";
void setVerboseLevel( const std::string& verboseLevel ) {
    m_debugLevel = verboseLevel;
    return;
}

void DEBUG(const std::string& message){
    /* Debug level (verbosity of output) */
    verbose("DEBUG",message);
    return;
}
void INFO(const std::string& message){
    /* Info level (verbosity of output) */
    verbose("INFO",message);
    return;
}
void WARNING(const std::string& message){
    /* Warning level (verbosity of output) */
    verbose("WARNING",message);
    return;
}
void ERROR(const std::string& message){
    /* Error level (verbosity of output) */
    verbose("ERROR",message);
    return;
}

void verbose(const std::string level, const std::string& message){
    /* 
       Printing output to console (debug,warning,error messages)
         if the level is "DEBUG", then all messages should be printed (DEBUG/INFO/WARNING/ERROR)
         if the level is "INFO", then only INFO/WARNING/ERROR messages should be printed
         if the level is "WARNING", then only WARNING/ERROR messages should be printed
         if the level is "ERROR", then only ERROR messages should be printed
    */
    std::map<std::string,unsigned int> debugMap = {
            {"DEBUG",   0},
            {"INFO",    1},
            {"WARNING", 2},
            {"ERROR",   3} };

    if ( debugMap.at( level ) >= debugMap.at( m_debugLevel ))
        std::cout << " " << level << " :: " << message << std::endl;

    return;
}

std::map<std::string,unsigned int> verboseMap() {
    /* mapping of verbose level to integer */
    std::map<std::string,unsigned int> verbose_map = {
            {"DEBUG",   0},
            {"INFO",    1},
            {"WARNING", 2},
            {"ERROR",   3} };
    
    return verbose_map;
}

void HELP(const std::string& runExecutable){
    /* HELP message (pass 'runExecutable' in case you are running from some 
       script like 'skim', 'run', or a custom macro)
    */
    std::cout << "\n   ** CyMiniAna ** " << std::endl;
    std::cout << "   --------------- " << std::endl;
    std::cout << "   Framework to perform event selection, write-out" << std::endl;
    std::cout << "   a few histograms or efficiencies, and make plots.\n" << std::endl;

    std::cout << "   To run:" << std::endl;
    std::cout << "      ./" << runExecutable << " share/cmaConfig.txt \n" << std::endl;
    std::cout << "    where 'share/cmaConfig.txt' is the configuration file \n" << std::endl;

    return;
}

} // end namespace

// the end
