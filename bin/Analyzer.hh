// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef ANALYZER_HH
#define ANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
//#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
//#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/tgPhysObject.hh"


// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class DemoAnalyzer: public BLTSelector {
public:
    DemoAnalyzer();
    ~DemoAnalyzer();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();

    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile       *outFile;
    TTree       *outTree;
		TTree 			*treeVector;
    std::string  outFileName;
    std::string  outTreeName;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
//    std::unique_ptr<Cuts>               cuts;
    std::unique_ptr<TriggerSelector>    triggerSelector;
//    std::unique_ptr<ParticleSelector>   particleSelector;

    //ClassDef(DemoAnalyzer,0);
};


#endif  // DEMOANALYZER_HH
