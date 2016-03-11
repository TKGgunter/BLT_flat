#include "DemoAnalyzer.hh"
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

int n_muons = 0;
int n_elec  = 0;
int gMuons = 0;
int gElectrons = 0;
int npv = 0;
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

//		std::cout << tmp_option << std::endl;
	
//    std::regex re_whitespace("(\\s+)");  // split by white space
//    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
//              std::sregex_token_iterator(), std::back_inserter(options));

    //std::cout << options.size() << std::endl;
    // Set the parameters
    //params.reset(new Parameters());
    //params->setup(options);

    // Set the cuts
//    cuts.reset(new Cuts());
//    particleSelector.reset(new ParticleSelector(*params, *cuts));
//    triggerSelector.reset(new TriggerSelector());

    // Prepare the output tree
    outFileName =  options.at(0) + "_" + options.at(2)  +"_" + options.at( options.size() - 1) + ".root";//"WJ_test.root";//params->get_output_filename("demoFile");
    outTreeName = "WJ";// params->get_output_treename("demoTree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
//    outTree = new TTree(outTreeName.c_str(), "demoTree");

//    outTree->Branch("muonOne", &muonOne);
//    outTree->Branch("muonTwo", &muonTwo);
//    outTree->Branch("dimuon", &dimuon);
//    outTree->Branch("genMuonOne", &genMuonOne);
//    outTree->Branch("genMuonTwo", &genMuonTwo);
//    outTree->Branch("genZ", &genZ);
//
  treeVector = new TTree("trees_vec", "simpTree_vec");

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
  treeVector->Branch("npv", &npv, "npv/I");


/**/
    ReportPostBegin();
}

Bool_t DemoAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
		//std::cout << "event # " << this->totalEvents << std::endl;
    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    //if (entry%1==0)  std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;

//    const bool isRealData = (fInfo->runNum != 1);
/*
    particleSelector->SetRealData(isRealData);
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        particleSelector->SetPV(TVector3());
    }
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);
    triggerSelector->SetRealData(isRealData);
*/
//    bool printEvent = false;
/*
		print_event(printEvent, isRealData, fMuonArr, fElectronArr, 
																				fTauArr, fAK4CHSArr,
																				fPhotonArr, fGenParticleArr,
																				fPVArr);
*/
/*
    //////////////////
    // GenParticles //
    //////////////////

    if (printEvent) {
        if (!isRealData) {
            for (int i=0; i<fGenParticleArr->GetEntries(); i++) {
                const TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
                assert(particle);
                std::cout << "GenParticle " << i << ": " << particle << std::endl;
            }
        }
    }

    //////////////
    // Vertices //
    //////////////

    if (printEvent) {
        for (int i=0; i<fPVArr->GetEntries(); i++) {
            const TVertex* vertex = (TVertex*) fPVArr->At(i);
            assert(vertex);
            std::cout << "Vertex " << i << ": " << vertex << std::endl;
        }
    }

    ///////////
    // Muons //
    ///////////

    if (printEvent) {
        for (int i=0; i<fMuonArr->GetEntries(); i++) {
            const TMuon* muon = (TMuon*) fMuonArr->At(i);
            assert(muon);
            std::cout << "Muon " << i << ": " << muon
                      << " tightMuID: " << particleSelector->PassMuonID(muon, cuts->tightMuID)
                      << " tightMuIso: " << particleSelector->PassMuonIso(muon, cuts->tightMuIso)
                      << std::endl;
        }
    }

    ///////////////
    // Electrons //
    ///////////////

    if (printEvent) {
        for (int i=0; i<fElectronArr->GetEntries(); i++) {
            const TElectron* electron = (TElectron*) fElectronArr->At(i);
            assert(electron);
            std::cout << "Electron " << i << ": " << electron << std::endl;
        }
    }

    /////////////
    // Photons //
    /////////////

    if (printEvent) {
        for (int i=0; i<fPhotonArr->GetEntries(); i++) {
            const TPhoton* photon = (TPhoton*) fPhotonArr->At(i);
            assert(photon);
            std::cout << "Photon " << i << ": " << photon << std::endl;
        }
    }

    //////////
    // Taus //
    //////////

    if (printEvent) {
        for (int i=0; i<fTauArr->GetEntries(); i++) {
            const TTau* tau = (TTau*) fTauArr->At(i);
            assert(tau);
            std::cout << "Tau " << i << ": " << tau << std::endl;
        }
    }

    //////////
    // Jets //
    //////////

    if (printEvent) {
        for (int i=0; i<fAK4CHSArr->GetEntries(); i++) {
            const TJet* jet = (TJet*) fAK4CHSArr->At(i);
            assert(jet);
            std::cout << "Jet " << i << ": " << jet << std::endl;
        }
    }
*/
    /////////
    // MET //
    /////////
 
	  TMET* pfMET = new TMET();
    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;

/*
    TMET* caloMET = new TMET();
    caloMET->pt = fInfo->caloMET;
    caloMET->phi = fInfo->caloMETphi;

    if (printEvent) {
        std::cout << "MET " << "(PF)" << ": " << pfMET << std::endl;
        std::cout << "MET " << "(C) " << ": " << caloMET << std::endl;
    }
*/
    ////////////
    // Select //
    ////////////


//////////////////////////////////////////////////////////////////////////////
/////////////////// Number of Lep Cut & vertex cut ///////////////////////////

		if(fElectronArr->GetEntries() + fMuonArr->GetEntries() < 2 ) return kTRUE; 


		//gMuons = 0;
		//gElectrons = 0;
		/*
		for(int i = 0; i < fGenParticleArr->GetEntries(); i++){
			TGenParticle* genPart = (TGenParticle*)fGenParticleArr->At(i);
			if ( abs(genPart->pdgId) == 11 && genPart->status == 1){
				//TGenParticle* motherPart = (TGenParticle*)fGenParticleArr->At(genPart->parent);
				//if(abs(motherPart->pdgId) == 24 )
				gElectrons++;
			}
			if ( abs(genPart->pdgId) == 13 &&  genPart->status == 1){
				//TGenParticle* motherPart = (TGenParticle*)fGenParticleArr->At(genPart->parent);
				//if(abs(motherPart->pdgId) == 24 )
				gMuons++;
			}
		}
	*/
	//&& eventConstrType == RECO) return kTRUE;

		int tallyVtx = 0;
		for(int i = 0; i < fPVArr->GetEntries() ; i++){
			TVertex* primVtx = (TVertex*) fPVArr->At(i);
			//if(	!primVtx->IsFake() && primVtx->ndof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ) tallyVtx++;
			if( primVtx->ndof > 4. && fabs(primVtx->z) <= 24. ) tallyVtx++;
		}

		if(tallyVtx == 0) return kTRUE;

////////////////////////////////////////////////
//////////////Muons Stuff//
		n_muons = n_muons + fMuonArr->GetEntries();
		n_elec = n_elec + fElectronArr->GetEntries();

		std::vector<TMuon*> muonList;
		//if(RECO)
		muonSelection(muonList, fMuonArr);

/*		for(int i = 0; i < fMuonArr->GetEntries(); i++){
			TMuon* muon = (TMuon*) fMuonArr->At(i);
			bool goodMu = particleSelector->PassMuonID(muon, cuts->tightMuID) && particleSelector->PassMuonIso(muon, cuts->tightMuIso);
			if( goodMu == true) muonList.push_back(muon);
		}
*/

		tgSort(muonList);

		if(muonList.size() != 0){
			tgCleanVector(muonList);
		}

/////////////Electron Stuff//
		std::vector<TElectron*> elecList;
		electronSelection(elecList, fElectronArr);

/*		for(int i = 0; i < fElectronArr->GetEntries(); i++){
			TElectron* electron = (TElectron*) fElectronArr->At(i);
			bool goodEle = particleSelector->PassElectronID(electron, cuts->looseElID) ;//&& particleSelector->PassElectronIso(electron, cuts->mediumElIso);
			if( goodEle == true) elecList.push_back(electron);
		}
*/

    tgSort(elecList);

    if(elecList.size() != 0){
      tgCleanVector(elecList);
    }

	//	int nLeps = elecList.size() + muonList.size();
		if(elecList.size() + muonList.size() < 2) return kTRUE;
//
//////////////Jet Stuff//
	  bool isJetGood = true;
  	std::vector<TJet*> jetList;

		for(Int_t i = 0; i < fAK4CHSArr->GetEntries(); i++){ 
//			TCJet* jet = (TCJet*) recoJets->At(i);
			TJet* jet = (TJet*) fAK4CHSArr->At(i);

//**********************//

//			TJet* testJet = jet_;
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
*/

      isJetGood = true;
      if( (jet->eta < jetEta) &&
          (jet->neuHadFrac < jetNeuHadFrac) &&
          (jet->neuEmFrac < jetNeuEmFrac) &&
          (jet->nParticles > jetNumConstit ) &&
          (jet->chEmFrac < jetChEmFrac) &&
          (jet->chHadFrac > jetChHadFrac) &&
          jet->pt > jetPt &&
          jet->nCharged > 0 )//jet->NumChPart() > 0 )
      {
 
      /* 
				if (jet->eta() < 2.5){
					if( jet->dR2Mean() > 0.025 || jet->BetaStarClassic() > 0.8 ) isJetGood = false;
        }
				if (jet->eta() > 2.75 && jet->eta() < 3 && jet->dR2Mean() > 0.04 ) isJetGood = false;
 
		*/

        if( elecList.size() == 0 ){
          for(  std::vector<TMuon*>::iterator muon = muonList.begin(); muon != muonList.end(); ++muon){
            float dR_muon = dR(*muon, jet);
            if(dR_muon < jetDR_me){
              isJetGood = false;
            }
          }
        }
        if(muonList.size() == 0){
          for(  std::vector<TElectron*>::iterator electron = elecList.begin(); electron != elecList.end(); ++electron){
            float dR_elec = dR(*electron, jet);
            if(dR_elec < jetDR_me){
              isJetGood = false;
            }
          }
        }
        if(muonList.size() != 0 && elecList.size() != 0){
          for(  std::vector<TMuon*>::iterator muon = muonList.begin(); muon != muonList.end(); ++muon){
            float dR_muon = dR(*muon, jet);
            if(dR_muon < jetDR_mu_el){
              isJetGood = false;
            }
          }
          for(  std::vector<TElectron*>::iterator electron = elecList.begin(); electron != elecList.end(); ++electron){
            float dR_elec = dR(*electron, jet);
            if(dR_elec < jetDR_mu_el){
              isJetGood = false;
            }
          }
        }
        if( isJetGood == true){
					jetList.push_back(jet);
				}
			}
		}
		tgSort(jetList);
		tgCleanVector(jetList , fPVArr->GetEntries());

	  //====================================================
	  //=================Set leaves ====================

	

		lep1_pt = 0.;
		lep2_pt = 0.;
		lep3_pt = 0.;
		lep_Charge = 0;	
		lep_Type = 0;
		dPhiLL = 0.;  
  	qT = 0.;
		dPhiLLJet = 0.;
		mll = 0.;
	
	//===============JETS=====================================
		jetPt1 = 0.;
		jetPt2 = 0.;
		jetPt3 = 0.;
		jetPt4 = 0.;
		jetPt5 = 0.;
		jetPt6 = 0.;
		numb_jets = 0;// numb_jets = jetList.size();


		numb_jets = jetList.size();
		if(jetList.size() > 0) jetPt1 = ((TJet*) jetList.at(0) )->pt;
		if(jetList.size() > 1) jetPt2 = ((TJet*) jetList.at(1) )->pt;
		if(jetList.size() > 2) jetPt3 = ((TJet*) jetList.at(2) )->pt;
		if(jetList.size() > 3) jetPt4 = ((TJet*) jetList.at(3) )->pt;
		if(jetList.size() > 4) jetPt5 = ((TJet*) jetList.at(4) )->pt;
		if(jetList.size() > 5) jetPt6 = ((TJet*) jetList.at(5) )->pt;
		
		nBJet = 0;
		HT = 0;
		for(unsigned int i = 0; i < jetList.size(); i++){
			TJet* jet = (TJet*) jetList.at(i);
			HT += jet->pt;
			//std::cout << jet->JetFlavor() << std::endl;
			if(abs(jet->mcFlavor) == 5) nBJet++; // Uses generator Flavor tag
		}

//======================Leps===========================
		std::vector<TGPhysObject*> particleList;
		tgConcateList(muonList, elecList, particleList);
		int flavor_mu = 0; int flavor_el = 0;

		for(unsigned int i = 0 ; i < particleList.size(); i++){
			TGPhysObject* particle = (TGPhysObject*) particleList.at(i);
			if(i == 0){	
				lep1_pt = particle->pt;
				lead_lep_type = particle->id;
				lead_lep_eta = particle->eta;
			}
			if(i == 1){	
				lep2_pt = particle->pt;
				sub_lep_type = particle->id;
				sub_lep_eta = particle->eta;
			}
			if(i == 2)	lep3_pt = particle->pt;

			if(i <= 1 && particle->id == 13) flavor_mu++;
			if(i <= 1 && particle->id == 11) flavor_el++;
		}

  	TGPhysObject* leadObj = (TGPhysObject*) particleList.at(0);
  	TGPhysObject* subObj = (TGPhysObject*) particleList.at(1);

  	qT = qtCalc(leadObj, subObj);
  	mll = ( (leadObj->p4) + (subObj->p4)).Mag();


  	dPhiLL = deltaPhi(subObj->phi - leadObj->phi);
  	dPhiLL = fabs(dPhiLL);
  	dPhiLLMET = dPhiLLJet = deltaPhi((subObj->p4 + leadObj->p4).Phi() - fInfo->pfMETphi);//recoMET->phi());
  	dPhiLLMET = fabs(dPhiLLMET);
  	if(jetList.size() != 0) dPhiLLJet = deltaPhi((subObj->p4 + leadObj->p4).Phi() - ((TJet*) jetList.at(0))->phi);
  	dPhiLLJet = fabs(dPhiLLJet);

  //lep_Charge = leadObj->Charge() + subObj->Charge();

  	lep_Type = 0;
  	if( flavor_mu == 2) lep_Type = -1 ;
  	if( flavor_el == 2) lep_Type = -2;
  	if( flavor_mu == 1 && flavor_el == 1 ) lep_Type = 1;

  	numbExtraLep = 0; numbExtraLep = muonList.size() + elecList.size() - 2;


//MET

    pfMET->pt = fInfo->pfMET;
    pfMET->phi = fInfo->pfMETphi;


	  float dPhi_m = 99;
  	float dPhi_e = 99;
  	TMuon* muon_minDr = NULL;
  	TElectron* electron_minDr = NULL;
  	int flagged_e = 0;
  	int flagged_m = 0;
  	for (unsigned int i = 0 ; i < muonList.size(); i++){
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

  	if(dPhi_m < dPhi_e ){
    	muon_minDr = (TMuon*) muonList.at(flagged_m);
    	if( dPhi_m < M_PI/2. ){
      	float muonMag = sqrt(pow(muon_minDr->pt * sin(muon_minDr->phi),2) + pow(muon_minDr->pt * cos(muon_minDr->phi),2) + pow(muon_minDr->pt * sinh(muon_minDr->eta),2));
      	METProj = sqrt(abs(pow(pfMET->pt, 2) - (pow(muon_minDr->pt*pfMET->pt * sin(muon_minDr->phi) * sin(pfMET->phi), 2) + pow(muon_minDr->pt * cos(muon_minDr->phi)*pfMET->pt * cos(pfMET->phi), 2))/pow(muonMag, 2)));
      	METProj_parr =  (pow(muon_minDr->pt*pfMET->pt * sin(muon_minDr->phi) * sin(pfMET->phi), 2) + pow(muon_minDr->pt*pfMET->pt * cos(muon_minDr->phi) * cos(pfMET->phi), 2))/pow(muonMag, 2);
    	}
    	else{
      	METProj = abs(pfMET->pt);
      	METProj_parr = abs(pfMET->pt);

    	}
  	}
  	else{
    	electron_minDr = (TElectron*) elecList.at(flagged_e);
    	if( dPhi_e < M_PI/2. ){
//      	float electronMag = sqrt(pow(electron_minDr->Px(),2) + pow(electron_minDr->Py(),2) + pow(electron_minDr->Pz(),2));         
//      	METProj = sqrt(abs(pow(pfMET->pt, 2) - (pow(electron_minDr->Px()*pfMET->Px(), 2) + pow(electron_minDr->Py()*pfMET->Py(), 2))/pow(electronMag, 2)));
//      	METProj_parr =  (pow(electron_minDr->Px()*pfMET->Px(), 2) + pow(electron_minDr->Py()*pfMET->Py(), 2))/pow(electronMag, 2);
        float electronMag = sqrt(pow(electron_minDr->pt * sin(electron_minDr->phi),2) + pow(electron_minDr->pt * cos(electron_minDr->phi),2) + pow(electron_minDr->pt * sinh(electron_minDr->eta),2));
        METProj = sqrt(abs(pow(pfMET->pt, 2) - (pow(electron_minDr->pt*pfMET->pt * sin(electron_minDr->phi) * sin(pfMET->phi), 2) + pow(electron_minDr->pt * cos(electron_minDr->phi)*pfMET->pt * cos(pfMET->phi), 2))/pow(electronMag, 2)));
        METProj_parr =  (pow(electron_minDr->pt*pfMET->pt * sin(electron_minDr->phi) * sin(pfMET->phi), 2) + pow(electron_minDr->pt*pfMET->pt * cos(electron_minDr->phi) * cos(pfMET->phi), 2))/pow(electronMag, 2);

      
    	}
    	else{
      	METProj = abs(pfMET->pt);
      	METProj_parr = abs(pfMET->pt);
    	}
  	}

  	metMod = 0.; metMod = pfMET->pt;

    //////////
    // Fill //
    //////////

    //muonOne = tmp_muonOne;

//    outTree->Fill();
		npv = tallyVtx;
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
