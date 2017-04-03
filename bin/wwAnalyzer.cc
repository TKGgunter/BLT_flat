#include "wwAnalyzer.hh"//#include "DemoAnalyzer.hh"
#include "BLT/BLTAnalysis/interface/SelectionCuts.hh"
#include "BLT/BLTAnalysis/interface/tgUtility.hh"
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
int lep1_Charge= 0;
int lep2_Charge= 0;
int lep_Type	= 0; 
int lep1_type = 0;
int lep2_type = 0;
float lep1_eta = 0;
float lep2_eta = 0;
float lep1_phi = 0;
float lep2_phi = 0;
float lep1_iso = 0;
float lep2_iso = 0;
int lep1_mother = 0;
int lep2_mother = 0;
float soft_muon_pt = 0;
float soft_muon_eta = 0;
float soft_muon_phi = 0;
float met_corrected	= 0.;
float met_trk	= 0.;
float metMod	= 0.;
float jetPt1	= 0.;
float jetPt2	= 0.;
float jetPt3	= 0.;
float jetPt4	= 0.;
float jetPt5	= 0.;
float jetPt6	= 0.;
float jet1_csv	= 0.;
float jet2_csv	= 0.;
float jet3_csv	= 0.;
float jet1_phi	= 0.;
float jet1_eta	= 0.;
float jet2_phi	= 0.;
float jet2_eta	= 0.;
int numbExtraLep =0;
float mll			= 0.;
float mllMET = 0.;
float recoil = 0.;
float dPhiLLJet = 0.;
int nBJet		= 0;
int nBJet_gen = 0;
float dPhiLLMET = 0.;
float dPhiMETJet = 0.;
float METProj 	= 0.;
float METProj_sin = 0.;
float METProj_trk_sin = 0.;
float MET_phi = 0.;
float met_corrected_phi = 0.;
float met_trk_phi = 0.;
float met_over_sET = 0.;
int nPartons = 0;
int tot_npv = 0;

int n_muons = 0;
int n_elec  = 0;
int gMuons = 0;
int gElectrons = 0;
int npv = 0;
float gen_npv = 0;

int runNum = 0;
int lumiSec = 0;
int eventNumb = 0;

float weight = 0.;
float npv_weight = 0;
float id_weight = 0;
std::string process; 
std::string process_decay; 

std::unique_ptr<WeightUtils>  weights;

using namespace baconhep;

DemoAnalyzer::DemoAnalyzer() : BLTSelector()
{

}

DemoAnalyzer::~DemoAnalyzer()
{

}



void DemoAnalyzer::Begin(TTree *tree)
{
	const std::string cmssw_base = getenv("CMSSW_BASE");
	weights.reset(new WeightUtils("a","a",false) );


	//***********************************//
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
	//**********************************//
	//**********************************//



	//*********************************//
	// Prepare the output tree
	outFileName =  options.at(0) + "_" + options.at(2)  +"_" + options.at( options.size() - 1) + ".root";//"WJ_test.root";//params->get_output_filename("demoFile");
	outTreeName = "WJ";// params->get_output_treename("demoTree");

	outFile = new TFile(outFileName.c_str(),"RECREATE");
	outFile->cd();
	
	printf("outfile set up\n");
	//********************************//
	//********************************//



	//********************************//
	// Set-up Triggers 
	triggerSelector.reset(new baconhep::TTrigger( cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2" ));
	// Set-up LumiMask
	lumiMask = RunLumiRangeMap();
	std::string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
	lumiMask.AddJSONFile( jsonFileName );

	printf("triggers set up\n");
	//******************************//
	//******************************//




	muonCorr = new rochcor2012();

	nEvents = new TH1F("nEvents", "nEvents", 10, 0, 10);
	nPV_plot = new TH1F("nPV", "nPV", 50, -0.5, 50.5);
	gen_nPV_plot = new TH1F("gen_nPV", "gen_nPV", 50, -0.5, 50.5);

  treeVector = new TTree("trees_vec", "simpTree_vec");

  treeVector->Branch("lep1_pt", &lep1_pt, "lep1_pt/F");
  treeVector->Branch("lep2_pt", &lep2_pt, "lep2_pt/F");
  treeVector->Branch("lep3_pt", &lep3_pt, "lep23_pt/F");
  treeVector->Branch("numb_jets", &numb_jets, "numb_jets/I");
  treeVector->Branch("dPhiLL", &dPhiLL, "dPhiLL/F");
  treeVector->Branch("qT", &qT, "qT/F");
  treeVector->Branch("HT", &HT, "HT/F");
  treeVector->Branch("lep1_Charge", &lep1_Charge,"lep1_Charge/I");
  treeVector->Branch("lep2_Charge", &lep2_Charge,"lep2_Charge/I");
  treeVector->Branch("lep_Type", &lep_Type, "lep_Type/I");
  treeVector->Branch("lep1_type", &lep1_type, "lep1_type/I");
  treeVector->Branch("lep2_type", &lep2_type, "lep2_type/I");
  treeVector->Branch("lep1_eta", &lep1_eta, "lep1_eta/F");
  treeVector->Branch("lep2_eta", &lep2_eta, "lep2_eta/F");
  treeVector->Branch("lep1_phi", &lep1_phi, "lep1_phi/F");
  treeVector->Branch("lep2_phi", &lep2_phi, "lep2_phi/F");
  treeVector->Branch("lep1_iso", &lep1_iso, "lep1_iso/F");
  treeVector->Branch("lep2_iso", &lep2_iso, "lep2_iso/F");
  treeVector->Branch("soft_muon_pt", &soft_muon_pt, "soft_muon_pt/F");
  treeVector->Branch("soft_muon_phi", &soft_muon_phi, "soft_muon_phi/F");
  treeVector->Branch("soft_muon_eta", &soft_muon_eta, "soft_muon_eta/F");
  treeVector->Branch("lep1_mother", &lep1_mother, "lep1_mother/I");
  treeVector->Branch("lep2_mother", &lep2_mother, "lep2_mother/I");
  treeVector->Branch("met_corrected", &met_corrected ,"met_corrected/F");
  treeVector->Branch("met_trk", &met_trk ,"met_trk/F");
  treeVector->Branch("metMod", &metMod ,"metMod/F");
  treeVector->Branch("jet1_pt", &jetPt1, "jet1_pt/F");
  treeVector->Branch("jet2_pt", &jetPt2, "jet2_pt/F");
  treeVector->Branch("jet3_pt", &jetPt3, "jet3_pt/F");
  treeVector->Branch("jet4_pt", &jetPt4, "jet4_pt/F");
  treeVector->Branch("jet5_pt", &jetPt5, "jet5_pt/F");
  treeVector->Branch("jet6_pt", &jetPt6, "jet6_pt/F");
  treeVector->Branch("jet1_csv", &jet1_csv, "jet1_csv/F");
  treeVector->Branch("jet2_csv", &jet2_csv, "jet2_csv/F");
  treeVector->Branch("jet3_csv", &jet3_csv, "jet3_csv/F");
  treeVector->Branch("jet1_phi", &jet1_phi, "jet1_phi/F");
  treeVector->Branch("jet2_phi", &jet2_phi, "jet2_phi/F");
  treeVector->Branch("jet1_eta", &jet1_eta, "jet1_eta/F");
  treeVector->Branch("jet2_eta", &jet2_eta, "jet2_eta/F");
  treeVector->Branch("numbExtraLep", &numbExtraLep, "numbExtraLep/I");
  treeVector->Branch("mll", &mll, "mll/F");
  treeVector->Branch("mllMET", &mllMET, "mllMET/F");
  treeVector->Branch("recoil", &recoil, "recoil/F");
  treeVector->Branch("dPhiLLJet", &dPhiLLJet, "dPhiLLJet/F");
  treeVector->Branch("numb_BJet", &nBJet, "numb_BJet/I");
  treeVector->Branch("numb_BJet_gen", &nBJet_gen, "numb_BJet_gen/I");
  treeVector->Branch("dPhiLLMET", &dPhiLLMET, "dPhiLLMET");
  treeVector->Branch("dPhiMETJet", &dPhiMETJet, "dPhiMETJet");
  treeVector->Branch("METProj", &METProj, "METProj");
  treeVector->Branch("METProj_sin", &METProj_sin, "METProj_sin");
  treeVector->Branch("METProj_trk_sin", &METProj_trk_sin, "METProj_trk_sin");
  treeVector->Branch("met_phi", &MET_phi, "met_phi");
  treeVector->Branch("met_corrected_phi", &met_corrected_phi, "met_corrected_phi");
  treeVector->Branch("met_trk_phi", &met_trk_phi, "met_trk_phi");
  treeVector->Branch("met_over_sET", &met_over_sET, "met_over_sET");
  treeVector->Branch("nPartons", &nPartons, "nPartons/I");

  treeVector->Branch("npv", &npv, "npv/I");
  treeVector->Branch("tot_npv", &tot_npv, "tot_npv/I");
  treeVector->Branch("gen_npv", &gen_npv, "gen_npv/F");
  treeVector->Branch("runNum", &runNum, "runNum/I");
  treeVector->Branch("lumiSec", &lumiSec, "lumiSec/I");
  treeVector->Branch("eventNumb", &eventNumb, "eventNumb/I");
	treeVector->Branch("weight", &weight, "weight/F" );
	treeVector->Branch("nov_weight", &npv_weight, "nov_weight/F" );
	treeVector->Branch("id_weight", &id_weight, "id_weight/F" );
	treeVector->Branch("process_decay", &process_decay);
	treeVector->Branch("process", &process);
	

	process_decay = options.at(0);	
	process = options.at(1) ;

	std::cout << " Process Decay stuff "  << process << " " << process_decay << std::endl;
/**/
    ReportPostBegin();

	printf("branches set up\n");

}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;



		int nEvents_bin = 0;
		nEvents->Fill( nEvents_bin ); nEvents_bin++;


		// SET UP TRIGGERS AND REAL DATA BOL
	//	printf("Triggers\n");
    const bool isRealData = (fInfo->runNum != 1);

		//Why is the following commented out
		//triggerSelector->SetRealData(isRealData);

		std::vector<std::string> trigger_vector;
		trigger_vector.push_back("HLT_IsoMu24_eta2p1_v*");

		trigger_vector.push_back("HLT_Ele27_WP80_v*");
//		trigger_vector.push_back("HLT_Mu17_v*");
		bool pass_trigger = false;
		for(unsigned i = 0; i < trigger_vector.size(); i++){
			pass_trigger |= triggerSelector->pass( trigger_vector[i], fInfo->triggerBits);
		}
		if( pass_trigger == false ) {
			//printf("trigger");
			return kTRUE;
		}


	//LUMI NUMBER CHECK
		//printf("LUMI\n");
		if (isRealData == 1){ 
			process = "Da";
			process_decay = "Da";
			RunLumiRangeMap::RunLumiPairType rl( fInfo->runNum, fInfo->lumiSec);
			if(!lumiMask.HasRunLumi(rl)) return kTRUE;
		}

		nEvents->Fill( nEvents_bin ); nEvents_bin++;


	//SET UP WEIGHTS
		//printf("weights start \n");
		weight = 1;
		npv_weight = 1;
		if (!isRealData){ 
			weight *= weights->GetPUWeight(fInfo->nPUmean);
			npv_weight *= weights->GetPUWeight(fInfo->nPUmean);
			//std::cout << "weights " << weight << std::endl;
		}


		//triggerSelector->SetTriggerNames(trigger_vector);

/*
		if( fInfo->runNum == 191226 && fInfo->lumiSec == 180 && fElectronArr->GetEntries() != 0){
			printf("Run Number 193621 lumiSec 804 Trigger test\n");
			printf("****************\n");
			for( int i = 0; i < fMuonArr->GetEntries(); i++){
				TMuon* muon = (TMuon*) fMuonArr->At(i);
				std::cout << "muon " << muon << " iso " << (muon->chHadIso03 + muon->neuHadIso03 + muon->gammaIso03) / muon->pfPt << " " << muon->typeBits << " " <<  kGlobal << " " << kPFMuon << std::endl;
			}
			for( int i = 0; i < fElectronArr->GetEntries(); i++){
				TElectron* elec = (TElectron*) fElectronArr->At(i);
				std::cout << "elec " << elec << " event " << fInfo->evtNum << std::endl;
			}
		}
		*/
		nEvents->Fill( nEvents_bin ); nEvents_bin++;



    /////////
    // MET //
    /////////
 
	  TMET* pfMET = new TMET();
    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;

	  TMET* trkMET = new TMET();
    trkMET->pt = fInfo->pfMET;
    trkMET->phi = fInfo->trkMETphi;

		TGPhysObject tgMET;
		tgMET.pt = pfMET->pt;
		tgMET.phi = pfMET->phi;
		tgMET.charge = 0;
		tgMET.id = 0;
		tgMET.calcP4();

    ////////////
    // Select //
    ////////////

//		if( fInfo->runNum == 191226 && fInfo->lumiSec == 180 && fElectronArr->GetEntries() != 0){
//			printf("TEST PASSED\n");
//		}
//////////////////////////////////////////////////////////////////////////////
/////////////////// Number of Lep Cut & vertex cut ///////////////////////////

		if(fElectronArr->GetEntries() + fMuonArr->GetEntries() < 2 ){
			return kTRUE; 
		}

		nEvents->Fill( nEvents_bin ); nEvents_bin++;
		
		int tallyVtx = 0;
		for(int i = 0; i < fPVArr->GetEntries() ; i++){
			TVertex* primVtx = (TVertex*) fPVArr->At(i);
			//if(	!primVtx->IsFake() && primVtx->ndof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ) tallyVtx++;
			if( primVtx->ndof > 4. && fabs(primVtx->z) <= 24. ) tallyVtx++;
		}


		
		if( fInfo->hasGoodPV == false ){ 
			return kTRUE;
		}
		nEvents->Fill( nEvents_bin ); nEvents_bin++;

////////////////////////////////////////////////
//////////////Muons Stuff//

		n_muons = n_muons + fMuonArr->GetEntries();
		n_elec = n_elec + fElectronArr->GetEntries();

		std::vector<TMuon*> muonList;
		std::vector<TMuon*> softMuonList;
		std::vector<TGenParticle*> muonFromTaus;
		if( isRealData == 1 ){ 
			muonSelection(muonList, fMuonArr);
			softMuonSelection(softMuonList, fMuonArr);
		}
		else{
			muonSelection(muonList, fMuonArr);
			softMuonSelection(softMuonList, fMuonArr);
			muonSelection(muonFromTaus, fGenParticleArr, 23, muonList);
		}
		tgSort(muonList);

		if(muonList.size() != 0){
			tgCleanVector(muonList);
		}
		if(softMuonList.size() != 0){
			tgCleanVector(softMuonList);
		}


/////////////Electron Stuff//

		std::vector<TElectron*> elecList;
		std::vector<TGenParticle*> elecFromTaus;
		if( isRealData == 1 ) electronSelection(elecList, fElectronArr, fInfo->rhoJet);
		else{
			electronSelection(elecList, fElectronArr, fInfo->rhoJet);
			electronSelection(elecFromTaus, fElectronArr, 23, elecList);
		}

    tgSort(elecList);

    if(elecList.size() != 0){
      tgCleanVector(elecList);
    }
//		if( fInfo->runNum == 191226 && fInfo->lumiSec == 180 ){
//			std::cout <<"!!!!!!!!!!!!muonlist size " << muonList.size() << " eleList size " << elecList.size() << " " << fInfo->evtNum << std::endl;
//		}

		if(elecList.size() + muonList.size() < 2){ 
			return kTRUE;
		}
		nEvents->Fill( nEvents_bin ); nEvents_bin++;

		if( fInfo->runNum == 191226 && fInfo->lumiSec == 180 && muonList.size() + elecList.size() == 3){
			printf("adfadfadf");
			for( unsigned int i = 0; i < muonList.size(); i++){
				TMuon* muon = (TMuon*) muonList.at(i);
				std::cout << "muon " << muon << " iso " << (muon->chHadIso03 + muon->neuHadIso03 + muon->gammaIso03) / muon->pfPt << " " << muon->trkIso03 << " " <<  kGlobal << " " << kPFMuon << std::endl;
			}
			for( unsigned int i = 0; i < elecList.size(); i++){
				TElectron* elec = (TElectron*) elecList.at(i);
				std::cout << "elec " << elec << " event " << fInfo->evtNum << std::endl;
			}
				
		}
/////////////////////////////
//			Check delta R between leptons
//
////////////////////////////




//////////////Jet Stuff//
	std::vector<TJet*> jetList;
	jetSelection( jetList, fAK4CHSArr, muonList, elecList, fPVArr->GetEntries());//Change to fAK5Arr maybe	
	tgSort(jetList);
	tgCleanVector(jetList);


//**********************//
/* Jet energy correction
 *
 *
			float sfEta[] = {1.052, 1.057, 1.096, 1.134, 1.288};
			float etaBins[] = {0., 0.5, 1.1, 1.7, 2.3, 5.0};	

			float matchPt = 0;
			int count = 0;	
			for(int j = 0; j < genJets->GetSize(); j++ ){
				TGenJet* genJet = (TCGenJet*) genJets->At(j);
				if(dR(testJet, genJet) < 0.15 ){
					matchPt += genJet->Pt();
					count++;
				}	
			}
			if(count > 0){
				for(int k = 0 ; k < 5; ++k){
					if( fabs(testJet->Eta()) > etaBins[k] && fabs(testJet->Eta()) < etaBins[k+1]){
						if( (matchPt -testJet->Pt()) / testJet->Pt() < 0.3 * testJet->Pt() ){
							float correctedPt = matchPt + sfEta[k] * (testJet->Pt() - matchPt );
							testJet->SetPtEtaPhiE( correctedPt, jet_->Eta(), jet_->Phi(), jet_->E()  );
						}
					}	
				}
			}
 
			TJet* jet = testJet;
*****************************/

	  //====================================================
	  //=================Set leaves ====================
	  //====================================================

		lep1_pt = 0.;
		lep2_pt = 0.;
		lep3_pt = 0.;
		lep1_Charge = 0;	
		lep2_Charge = 0;	
		lep_Type = 0;
		dPhiLL = 0.;  
  	qT = 0.;
		dPhiLLJet = 0.;
		dPhiMETJet = 0.;
		mll = 0.;
		mllMET = 0.;
		recoil = 0.;
		soft_muon_pt = 0.;
		soft_muon_phi = 0.;
		soft_muon_eta = 0.;
	
	//===============JETS=====================================
		jetPt1 = 0.;
		jetPt2 = 0.;
		jetPt3 = 0.;
		jetPt4 = 0.;
		jetPt5 = 0.;
		jetPt6 = 0.;
		jet1_csv = -2.;
		jet2_csv = -2.;
		jet3_csv = -2.;
		numb_jets = 0;// numb_jets = jetList.size();


		numb_jets = jetList.size();
		if(jetList.size() > 0){
			jetPt1 = ((TJet*) jetList.at(0) )->pt;
			jet1_csv = ((TJet*) jetList.at(0) )->csv;
			jet1_phi = ((TJet*) jetList.at(0) )->phi;
			jet1_eta = ((TJet*) jetList.at(0) )->eta;
			dPhiMETJet = fabs( deltaPhi( fInfo->pfMETphi - ((TJet*) jetList.at(0) )->phi) );	
		}
		if(jetList.size() > 1){
			jetPt2 = ((TJet*) jetList.at(1) )->pt;
			jet2_csv = ((TJet*) jetList.at(1) )->csv;
			jet2_phi = ((TJet*) jetList.at(1) )->phi;
			jet2_eta = ((TJet*) jetList.at(1) )->eta;
		}
		if(jetList.size() > 2){
			jetPt3 = ((TJet*) jetList.at(2) )->pt;
			jet3_csv = ((TJet*) jetList.at(2) )->csv;
		}
		if(jetList.size() > 3) jetPt4 = ((TJet*) jetList.at(3) )->pt;
		if(jetList.size() > 4) jetPt5 = ((TJet*) jetList.at(4) )->pt;
		if(jetList.size() > 5) jetPt6 = ((TJet*) jetList.at(5) )->pt;

	
		nBJet = 0;
		nBJet_gen = 0;
		HT = 0;
		for(unsigned int i = 0; i < jetList.size(); i++){
			TJet* jet = (TJet*) jetList.at(i);
			HT += jet->pt;
			//std::cout << jet->JetFlavor() << std::endl;
			if(abs(jet->mcFlavor) == 5) nBJet_gen++; // Uses generator Flavor tag
			if( jet->csv > 0.898 ) nBJet++;
		}

//======================Leps===========================

		std::vector<TGPhysObject*> particleList;
		std::vector<TGPhysObject*> fromTausList;
		tgConcateList(muonList, elecList, particleList);
		tgConcateList(muonFromTaus, elecFromTaus, fromTausList);
	
		//muon energy correction	
		for(unsigned int i = 0; i < particleList.size(); i++ ){
			TGPhysObject* particle = (TGPhysObject*) particleList.at(i);
			if ( particle->id == 13){
				float qter = 1.;
				//std::cout << " before "  << particle->pt << " " << particle->phi << " " << particle->eta << std::endl;
				if (isRealData){
					muonCorr->momcor_data(particle->p4, float(particle->charge), 0, qter);
				}
				else muonCorr->momcor_mc(particle->p4, float(particle->charge), 0, qter);
			
				particle->decontructP4();
				//std::cout << " after "  << particle->pt << " " << particle->phi << " " << particle->eta << std::endl;
			}
			
		}	

		std::cout<<"soft muon "<<softMuonList.size()<< " reg muon " << muonList.size()<<std::endl;
		if ( softMuonList.size() > 0){
			TMuon* soft_muon = softMuonList.at(0);
			soft_muon_pt = soft_muon->pfPt;
			soft_muon_eta = soft_muon->pfEta;
			soft_muon_phi = soft_muon->pfPhi;
		//	float reg_muon_pt = ((TMuon*)muonList.at(0))->pfPt;
		//	std::cout << soft_muon_pt << " " << reg_muon_pt << " " << soft_muon_phi << std::endl;
		}

		if( isRealData != 1 ){ 
			if( fromTausList.size() >= 2 ){
				for(unsigned int i = 0; i < particleList.size(); i++ ){
					TGPhysObject* particle = (TGPhysObject*) particleList.at(i);

					for(unsigned int j = 0; j < fromTausList.size(); j++){
						TGPhysObject* fTaus = (TGPhysObject*) fromTausList.at(j);

						if( dR( particle, fTaus) < 0.3 ){ 
							particle->mother = 15;
						}
					}
				}
			}
		}

		int flavor_mu = 0; int flavor_el = 0;

		for(unsigned int i = 0 ; i < particleList.size(); i++){
			TGPhysObject* particle = (TGPhysObject*) particleList.at(i);
			if(i == 0){	
				lep1_pt = particle->pt;
				lep1_type = particle->id;
				lep1_eta = particle->eta;
				lep1_phi = particle->phi;
				lep1_mother = particle->mother;
				lep1_iso = particle->iso;
			}
			if(i == 1){	
				lep2_pt = particle->pt;
				lep2_type = particle->id;
				lep2_eta = particle->eta;
				lep2_phi = particle->phi;
				lep2_mother = particle->mother;
				lep2_iso = particle->iso;
			}
			if(i == 2)	lep3_pt = particle->pt;

			if(i <= 1 && particle->id == 13) flavor_mu++;
			if(i <= 1 && particle->id == 11) flavor_el++;
		}
		if (lep1_pt < 24) return kTRUE;

  	TGPhysObject* leadObj = (TGPhysObject*) particleList.at(0);
  	TGPhysObject* subObj = (TGPhysObject*) particleList.at(1);

		//Muon Weights
		id_weight = 1;
		if (!isRealData) {
			if(leadObj->id == 13) id_weight *= weights->GetMuonRecoEff(leadObj->p4);
			else id_weight *= weights->GetElectronRecoEff(leadObj->p4);
			if(subObj->id == 13) id_weight *= weights->GetMuonRecoEff(subObj->p4);
			else id_weight *= weights->GetElectronRecoEff(subObj->p4);
			weight *= id_weight;
		}

  	qT = qtCalc(leadObj, subObj);
  	mll = ( (leadObj->p4) + (subObj->p4)).Mag();
		mllMET = ( (leadObj->p4) + (subObj->p4) + tgMET.p4 ).Mag();
		recoil = ( (leadObj->p4).Vect() + (subObj->p4).Vect() + (tgMET.p4).Vect() ).Pt();
		//std::cout << "recoil " << recoil << " \n Px " << leadObj->p4.Px() << " " << subObj->p4.Px() << " " << tgMET.p4.Px() << std::endl;
		//std::cout << " Py " << leadObj->p4.Py() << " " << subObj->p4.Py() << " " << tgMET.p4.Py() << std::endl;
		//std::cout << " Pz " << leadObj->p4.Pz() << " " << subObj->p4.Pz() << " " << tgMET.p4.Pz() << std::endl;


  	dPhiLL = deltaPhi(subObj->phi - leadObj->phi);
  	dPhiLL = fabs(dPhiLL);
  	dPhiLLMET = deltaPhi((subObj->p4 + leadObj->p4).Phi() - fInfo->pfMETphi);//recoMET->phi());
  	dPhiLLMET = fabs(dPhiLLMET);
  	if(jetList.size() != 0) dPhiLLJet = deltaPhi((subObj->p4 + leadObj->p4).Phi() - ((TJet*) jetList.at(0))->phi);
  	dPhiLLJet = fabs(dPhiLLJet);

		lep1_Charge = leadObj->charge;
		lep2_Charge = subObj->charge;

  	lep_Type = 0;
  	if( flavor_mu == 2) lep_Type = -1 ;
  	if( flavor_el == 2) lep_Type = -2;
  	if( flavor_mu == 1 && flavor_el == 1 ) lep_Type = 1;

  	numbExtraLep = 0; 
		numbExtraLep = muonList.size() + elecList.size() - 2;


//MET

    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;


	  float dPhi_1 =  fabs( deltaPhi( leadObj->phi - pfMET->phi ) );
  	float dPhi_2 =  fabs( deltaPhi( subObj->phi - pfMET->phi ) );
	  float dPhi_trk_1 =  fabs( deltaPhi( leadObj->phi - trkMET->phi ) );
  	float dPhi_trk_2 =  fabs( deltaPhi( subObj->phi - trkMET->phi ) );
/*
 * for (unsigned int i = 0 ; i < muonList.size(); i++){
    	TMuon* muon = (TMuon*) muonList.at(i);
    	if(deltaPhi(muon->phi - pfMET->phi) < dPhi_m){
      	dPhi_m = deltaPhi(muon->phi-pfMET->phi);
      	flagged_m = i;
    	}
  	}
  	for (unsigned int i = 0 ; i < elecList.size() ; i++){
    	TElectron* electron = (TElectron*) elecList.at(i);
    	if(deltaPhi(electron->phi - pfMET->phi) < dPhi_e) {
      	dPhi_e = deltaPhi(electron->phi - pfMET->phi);
      	flagged_e = i;
    	}
  	}
*/

			METProj = 0.;
			METProj_sin = 0.;


			METProj_trk_sin = 0.;
			if( dPhi_trk_1 < dPhi_trk_2){
				if( abs(dPhi_trk_1) < M_PI/2. ){
					METProj_trk_sin =  fabs(sin(dPhi_trk_1)*trkMET->pt);
				}
				else METProj_trk_sin = trkMET->pt;
			}
			else{
				if( abs(dPhi_trk_2) < M_PI/2. ){
					METProj_trk_sin =  fabs(sin(dPhi_trk_2)*trkMET->pt);
				}
				else METProj_trk_sin = trkMET->pt;
			}


			if( dPhi_1 < dPhi_2 ){
				if( abs(dPhi_1) < M_PI/2. ){

					float muonMag = sqrt(pow(leadObj->pt * sin(leadObj->phi),2) + pow(leadObj->pt * cos(leadObj->phi),2) + pow(leadObj->pt * sinh(leadObj->eta),2));
					METProj = sqrt(abs(pow(pfMET->pt, 2) - (pow(leadObj->pt*pfMET->pt * sin(leadObj->phi) * sin(pfMET->phi), 2) + pow(leadObj->pt * cos(leadObj->phi)*pfMET->pt * cos(pfMET->phi), 2))/pow(muonMag, 2)));

					METProj_sin =  fabs(sin(dPhi_1)*pfMET->pt);
				}
				else{
					METProj = fabs(pfMET->pt);
					METProj_sin = fabs(pfMET->pt);

				}
			}
  	else{
    	if( abs(dPhi_2) < M_PI/2. ){

        float electronMag = sqrt(pow(subObj->pt * sin(subObj->phi),2) + pow(subObj->pt * cos(subObj->phi),2) + pow(subObj->pt * sinh(subObj->eta),2));
        METProj = sqrt(abs(pow(pfMET->pt, 2) - (pow(subObj->pt*pfMET->pt * sin(subObj->phi) * sin(pfMET->phi), 2) + pow(subObj->pt * cos(subObj->phi)*pfMET->pt * cos(pfMET->phi), 2))/pow(electronMag, 2)));

        METProj_sin =  fabs(sin(dPhi_2)*pfMET->pt);
    	}
    	else{
      	METProj = fabs(pfMET->pt);
      	METProj_sin = fabs(pfMET->pt);
    	}
  	}
		//if(abs(dPhi_m) < M_PI / 2.) std::cout << "orig " << METProj << "\tcos " << abs(sin(dPhi_m)*pfMET->pt) << "\tmet "<< pfMET->pt << "\tDphi " <<  dPhi_m << "" << std::endl;

  	metMod = 0.; metMod = pfMET->pt;
		met_trk = 0.; met_trk = trkMET->pt; 
		met_corrected = 0.; met_corrected = fInfo->pfMETC; 
		met_corrected_phi = 0.; met_corrected_phi = fInfo->pfMETCphi; 
		met_trk_phi = 0.; met_trk_phi = trkMET->phi; 
		met_over_sET = 0.; met_over_sET = pfMET->pt / sqrt( pow((leadObj->p4).Et(),2.) + pow((subObj->p4).Et(),2.) ) ; 

    //////////
    // Fill //
    //////////

		nPartons = -999;
		if( !isRealData ){
			unsigned count = 0;
			for( int i = 0; i < fGenParticleArr->GetEntries(); i++){
				TGenParticle* particle = (TGenParticle*)fGenParticleArr->At(i);

				if( particle->status == 3 && (abs(particle->pdgId) < 6 || particle->pdgId == 21) ) ++count; 
			}
			nPartons = count-4;
		}
		else{
			nPartons = 0;
		}
		tot_npv = -999;
		tot_npv = fPVArr->GetEntries();

		MET_phi = pfMET->phi;
		npv = fInfo->nPU + 1;
		gen_npv = fInfo->nPUmean;
		runNum = fInfo->runNum;
		lumiSec = fInfo->lumiSec;
		eventNumb = fInfo->evtNum; 
		nEvents->Fill( nEvents_bin );
		nPV_plot->Fill( tot_npv, weight);
		gen_nPV_plot->Fill( gen_npv);

		treeVector->Fill();
    this->passedEvents++;

    delete pfMET;
  //  delete caloMET;

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
//    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
//    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
		std::cout << " RECO      : Muons "  << n_muons << "\t" << "Electrons " << n_elec << std::endl;
//		std::cout << " GEN       : Muons "  << gMuons << "\t" << "Electrons " << gElectrons << std::endl;
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
