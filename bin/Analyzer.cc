#include "DemoAnalyzer.hh"
//#include "BLT/BLTAnalysis/interface/SelectionCuts.hh"
//#include "BLT/BLTAnalysis/interface/tgUtility.hh"
#include <cmath>
// See header file for class documentation
//

float lep1_pt = 0.;
float lep2_pt = 0.;
float lep3_pt = 0.;
int numb_jets = 0;
float dPhiLL	= 0.;
float qT 			= 0.;
float HT			= 0.;
int lep_Charge= 0;
int lep_Type	= 0; 
int lead_lep_type = 0;
int sub_lep_type = 0;
float lead_lep_eta = 0;
float sub_lep_eta = 0;
float metMod	= 0.;
float jetPt1	= 0.;
float jetPt2	= 0.;
float jetPt3	= 0.;
float jetPt4	= 0.;
float jetPt5	= 0.;
float jetPt6	= 0.;
int numbExtraLep =0;
float mll			= 0.;
float dPhiLLJet = 0.;
int nBJet		= 0;
float dPhiLLMET = 0.;
float METProj 	= 0.;
float METProj_parr = 0.;
float MET_phi = 0.;

int n_muons = 0;
int n_elec  = 0;
int gMuons = 0;
int gElectrons = 0;
int npv = 0;

int runNum = 0;
int lumiSec = 0;
int eventNumb = 0;
/**/
using namespace baconhep;

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}

void DemoAnalyzer::Begin(TTree *tree)
{

    // Parse command line option
	std::string tmp_option = GetOption();
  std::vector<std::string> options;


	bool more_to_parse = true;
	int i = 0;
	int iter = 0;
	std::cout << tmp_option << std::endl;
	while(more_to_parse == true){
		std::string token = tmp_option.substr(i, tmp_option.find(" ", i, 1) - i );
		//std::cout << token << " " << i << " " << tmp_option.find(" ", i, 1) << std::endl;
		options.push_back(token);
		iter += 1;
		if( tmp_option.find(" ", i, 1) == std::string::npos ) more_to_parse = false;
		i = tmp_option.find(" ", i, 1) + 1 ;
	}


	// Prepare the output tree
	outFileName =  options.at(0) + "_" + options.at(2)  +"_NPV_" + options.at( options.size() - 1) + ".root";

	outFile = new TFile(outFileName.c_str(),"RECREATE");
	outFile->cd();

  treeVector = new TTree("trees_vec", "simpTree_vec");
/*
  treeVector->Branch("lep1_pt", &lep1_pt, "lep1_pt/F");
  treeVector->Branch("lep2_pt", &lep2_pt, "lep2_pt/F");
  treeVector->Branch("lep3_pt", &lep3_pt, "lep23_pt/F");
  treeVector->Branch("numb_jets", &numb_jets, "numb_jets/I");
  treeVector->Branch("dPhiLL", &dPhiLL, "dPhiLL/F");
  treeVector->Branch("qT", &qT, "qT/F");
  treeVector->Branch("HT", &HT, "HT/F");
  treeVector->Branch("lep_Charge", &lep_Charge,"lep_Charge/I");
  treeVector->Branch("lep_Type", &lep_Type, "lep_Type/I");
  treeVector->Branch("lead_lep_type", &lead_lep_type, "lead_lep_type/I");
  treeVector->Branch("sub_lep_type", &sub_lep_type, "sub_lep_type/I");
  treeVector->Branch("lead_lep_eta", &lead_lep_eta, "lead_lep_eta/F");
  treeVector->Branch("sub_lep_eta", &sub_lep_eta, "sub_lep_eta/F");
  treeVector->Branch("metMod", &metMod ,"metMod/F");
  treeVector->Branch("jetPt1", &jetPt1, "jetPt1/F");
  treeVector->Branch("jetPt2", &jetPt2, "jetPt2/F");
  treeVector->Branch("jetPt3", &jetPt3, "jetPt3/F");
  treeVector->Branch("jetPt4", &jetPt4, "jetPt4/F");
  treeVector->Branch("jetPt5", &jetPt5, "jetPt5/F");
  treeVector->Branch("jetPt6", &jetPt6, "jetPt6/F");
  treeVector->Branch("numbExtraLep", &numbExtraLep, "numbExtraLep/I");
  treeVector->Branch("mll", &mll, "mll/F");
  treeVector->Branch("dPhiLLJet", &dPhiLLJet, "dPhiLLJet/F");
  treeVector->Branch("nBJet", &nBJet, "nBJet/I");
  treeVector->Branch("dPhiLLMET", &dPhiLLMET, "dPhiLLMET");
  treeVector->Branch("METProj", &METProj, "METProj");
  treeVector->Branch("METProj_parr", &METProj_parr, "METProj_parr");
  treeVector->Branch("MET_phi", &MET_phi, "MET_phi");
*/
  treeVector->Branch("npv", &npv, "npv/I");
/*
  treeVector->Branch("runNum", &runNum, "runNum/I");
  treeVector->Branch("lumiSec", &lumiSec, "lumiSec/I");
  treeVector->Branch("eventNumb", &eventNumb, "eventNumb/I");
*/

/**/
    ReportPostBegin();
}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
		npv = fInfo->nPU + 1;
		treeVector->Fill();
    return kTRUE;
}

void DemoAnalyzer::Terminate()
{

    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void DemoAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
//    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void DemoAnalyzer::ReportPostTerminate()
{

    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  ============================================================" << std::endl;

}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<DemoAnalyzer> selector(new DemoAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
