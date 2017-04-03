#include "BLT/BLTAnalysis/interface/SelectionCuts.hh"
//#include "tgUtility.C"

bool test_bits_(unsigned int bits, unsigned int test) {
    return (bits & test) == test;
}
void jetSelection( std::vector<baconhep::TJet*> &jetList , TClonesArray* jets, std::vector<baconhep::TMuon*> &muonList, std::vector<baconhep::TElectron*> &elecList, int nVtx){

	
	for( Int_t i = 0; i < jets->GetEntries(); i++){
		baconhep::TJet* jet = (baconhep::TJet*) jets->At(i);
		bool isJetGood = true;
		if( (fabs(jet->eta) < jetEta) &&
				(jet->pt > jetPt) &&
				(jet->neuHadFrac < jetNeuHadFrac) &&
				(jet->neuEmFrac < jetNeuEmFrac) &&
				(jet->nParticles > jetNumConstit ) )
		{
			if( (fabs(jet->eta) < 2.4) &&
					!(jet->chEmFrac < jetChEmFrac) &&
					!(jet->nCharged > 0 ) && 
					!(jet->chHadFrac > jetChHadFrac) )  isJetGood = false;

			//deltaR cuts
			if( elecList.size() == 0 ){
				for(  std::vector<baconhep::TMuon*>::iterator muon = muonList.begin(); muon != muonList.end(); ++muon){
					float dR_muon = dR(*muon, jet);
					if(dR_muon < jetDR_me){
						isJetGood = false;
					}
				}
			}
			if(muonList.size() == 0){
				for(  std::vector<baconhep::TElectron*>::iterator electron = elecList.begin(); electron != elecList.end(); ++electron){
					float dR_elec = dR(*electron, jet);
					if(dR_elec < jetDR_me){
						isJetGood = false;
					}
				}
			}
			if(muonList.size() != 0 && elecList.size() != 0){
				for(  std::vector<baconhep::TMuon*>::iterator muon = muonList.begin(); muon != muonList.end(); ++muon){
					float dR_muon = dR(*muon, jet);
					if(dR_muon < jetDR_mu_el){
						isJetGood = false;
					}
				}
				for(  std::vector<baconhep::TElectron*>::iterator electron = elecList.begin(); electron != elecList.end(); ++electron){
					float dR_elec = dR(*electron, jet);
					if(dR_elec < jetDR_mu_el){
						isJetGood = false;
					}
				}
			}

			if(fabs(jet->eta) < 2.5){
				if ( jet->betaStar/log(nVtx-0.64) > .2  && jet->dR2Mean > 0.06) {
					isJetGood = false;
				}
			}
			else if(fabs(jet->eta) < 2.75  ){
				if (jet->betaStar/log(nVtx-0.64) > .3  && jet->dR2Mean > 0.05){ 
					isJetGood = false;
				}
			}
			else if( fabs( jet->eta) < 3){
				if ( jet->dR2Mean > 0.05) {
					isJetGood = false; 
				}
			}

			else {
				if( jet->dR2Mean > 0.055){
					isJetGood = false;
				}
			}


			if( isJetGood == true){
				jetList.push_back(jet);
			}


		}
	}
}

//Muon Selection
void muonSelection(std::vector<baconhep::TMuon*> &muonList, TClonesArray* muons){
	std::vector<baconhep::TMuon*> tmp_muons;
	for(Int_t i = 0; i < muons->GetEntries(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muons->At(i);
		float sumPfIso = (muon->chHadIso03 + muon->neuHadIso03 + muon->gammaIso03) / muon->pfPt;

		if (
				fabs(muon->pfEta) < 2.1 &&//muEta   &&
				//muon->IsGood() == muIsGood &&
				test_bits_(muon->typeBits, baconhep::kGlobal) == muIsGLB && 
				//test_bits_(muon->typeBits, baconhep::kPFMuon) == muIsPF &&
				int(muon->nMatchStn) > 1 && // muNumbMatchedStat &&
				int(muon->nTkLayers)  > 5 &&//muTrackLayersMeasurements &&
				muon->muNchi2 < 10 &&//muNormChi2 &&
				fabs(muon->d0) < 0.2 && //muDxy &&
				fabs(muon->dz) < 0.5 &&//muDz &&
				muon->nValidHits > 0 && // muNumbValidHits &&
				sumPfIso < muSumPFIso && //pfIso
				muon->pfPt > 15 //muPt
			)
		{
			
			//muonList.push_back(muon);

			//Track Iso
			tmp_muons.push_back(muon);

		}
	}
////////////////
	for( unsigned int i = 0; i < tmp_muons.size(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) tmp_muons.at(i);

		//std::cout << tmp_muons.size() << " i " << i << std::endl;
		for( unsigned int j = i+1; j < tmp_muons.size(); j++){
			baconhep::TMuon* muon_j = (baconhep::TMuon*) tmp_muons.at(j);

			if( dR( muon_j->phi, muon->phi, muon_j->eta, muon->eta) ){
				muon->trkIso03 = std::max((float)0., muon->trkIso03 - muon_j->pt);
				muon_j->trkIso03 = std::max((float)0., muon_j->trkIso03 - muon->pt);
			}
		}

		if( muon->trkIso03/muon->pt  < 0.1){ 
			muonList.push_back(muon); 
			//std::cout << "asdf " << muon->trkIso03  << " " << muon->trkIso03/muon->pt<< std::endl;
		}
	}
///////////////
} 

void softMuonSelection(std::vector<baconhep::TMuon*> &muonList, TClonesArray* muons){
	std::vector<baconhep::TMuon*> tmp_muons;
	for(Int_t i = 0; i < muons->GetEntries(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muons->At(i);
		//float sumPfIso = (muon->chHadIso03 + muon->neuHadIso03 + muon->gammaIso03) / muon->pfPt;

		if (
				fabs(muon->pfEta) < 2.4 &&//muEta   &&
				//muon->IsGood() == muIsGood &&
				//test_bits_(muon->typeBits, baconhep::kGlobal) == muIsGLB && 
				//test_bits_(muon->typeBits, baconhep::kPFMuon) == muIsPF &&
				int(muon->nMatchStn) >= 1 && // muNumbMatchedStat &&
				muon->muNchi2 < 1.8 &&//muNormChi2 &&
				fabs(muon->d0) < 3 && //muDxy &&
				fabs(muon->dz) < 30 &&//muDz &&
				muon->nValidHits > 10 && // muNumbValidHits &&
				muon->pfPt > 5 //muPt
			)
		{
			//Track Iso
			tmp_muons.push_back(muon);

		}
	}
////////////////
	for( unsigned int i = 0; i < tmp_muons.size(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) tmp_muons.at(i);

		//std::cout << tmp_muons.size() << " i " << i << std::endl;
		for( unsigned int j = i+1; j < tmp_muons.size(); j++){
			baconhep::TMuon* muon_j = (baconhep::TMuon*) tmp_muons.at(j);

			if( dR( muon_j->phi, muon->phi, muon_j->eta, muon->eta) ){
				muon->trkIso03 = std::max((float)0., muon->trkIso03 - muon_j->pt);
				muon_j->trkIso03 = std::max((float)0., muon_j->trkIso03 - muon->pt);
			}
		}

		if( muon->trkIso03/muon->pt  < 0.1){ 
			muonList.push_back(muon); 
			//std::cout << "asdf " << muon->trkIso03  << " " << muon->trkIso03/muon->pt<< std::endl;
		}
	}
///////////////
} 

void muonSelection(std::vector<baconhep::TGenParticle*> &muonList, TClonesArray* muons, int id, std::vector<baconhep::TMuon*> &recoMuonList ){

	if (id == 0) return;
	bool push_back = false;	
	bool wasTau = false;
	//printf("test1 muon selection\n genparticle size %i\n", muons->GetSize());
 
  for(Int_t i = 0; i < muons->GetSize(); i++){
    baconhep::TGenParticle* muon = (baconhep::TGenParticle*) muons->At(i);
		if( muon == NULL ) break;

		wasTau = false;
		//std::cout << "test2 " <<  muon->pdgId << " iter " << i << " " << push_back << wasTau <<std::endl;	

		if ( abs(muon->pdgId) == 13 && muon->status == 1){
			//printf("test3 parent %i \n", muon->parent);
			if(muon->parent < 0 ) break;
			if(muons->At( muon->parent ) == NULL) break;//return;
			baconhep::TGenParticle* mother;

			//std::cout << "muon parent " << muon->parent << std::endl;

			mother  = (baconhep::TGenParticle*)( muons->At( muon->parent ) );
			baconhep::TGenParticle* a;
			

			while(abs(mother->pdgId) == 13 || abs(mother->pdgId) == 11 || abs(mother->pdgId) == 15){
				//std::cout << "test4 mother id " <<  mother->pdgId  << " mother's parent " << mother->parent << std::endl;
				if( mother->parent < 0) break;
				if(muons->At( mother->parent ) == NULL) break;//return;
				//std::cout << "test4 mother id " <<  mother->pdgId  << " mother's parent " << mother->parent << std::endl;
				if( abs(mother->pdgId) == 15 ) wasTau = true;
				a = (baconhep::TGenParticle*)( muons->At( mother->parent ) );
				mother = a;
				
			}
			
			//if (abs(mother->GetPDGId()) != id) break;// != id return; 

			if ( abs(mother->pdgId) == id &&
				  wasTau == true						&&	
					fabs(muon->eta) < muEta   &&
		      //fabs(muon->Dz((baconhep::TVertex*)vtxList.at(0)))  < muDz &&
	        //fabs(muon->Dxy((baconhep::TVertex*)vtxList.at(0))) < muDxy &&
					muon->pt > muPt)
			{
				push_back = false;
				for( unsigned int j = 0; j < recoMuonList.size(); j++ ){
//					std::cout << j << std::endl;
					baconhep::TMuon* recoMuon = (baconhep::TMuon*) recoMuonList.at(j);
					if( dR( muon->phi , recoMuon->phi, muon->eta, recoMuon->eta) < 0.3 ){
						push_back = true;
					}
				}


				if( push_back == true){ 
					muonList.push_back((baconhep::TGenParticle*)muon);
					//printf("muon %f\n", muon->pt);
				}
			}
		}
  }

}



//Electron Selection
void electronSelection(std::vector<baconhep::TElectron*> &elecList, TClonesArray* electrons, float rhoFactor){

	int iEta = 0;
	float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
	float EAEL[7] = {0.13, 0.14, 0.07, 0.09, 0.11, 0.11, 0.14};

	for(Int_t i = 0; i < electrons->GetEntries(); i++){
		baconhep::TElectron* electron = (baconhep::TElectron*) electrons->At(i);
		float sumPfIso = (electron->chHadIso03 + electron->neuHadIso03 + electron->gammaIso03 )/electron->pfPt;
		float energyInversion = fabs( 1. - electron->eoverp) / electron->ecalEnergy;		

		for( unsigned int j = 0; j < 8; ++j ){
			if( fabs(electron->scEta) > etaBins[j] && fabs(electron->scEta) < etaBins[j+1] ){
				iEta = j;
			}
		}


		float combIso = electron->chHadIso03 + std::max(0., (double)electron->neuHadIso03 + electron->gammaIso03 - rhoFactor * EAEL[iEta]);
		sumPfIso = combIso / electron->pt;
	
		if(fabs(electron->scEta) < 1.479 ){//electron->IsEB() == true){

			if(	fabs(electron->eta) < elEta &&
					electron->pt > elPt &&
					energyInversion < 0.05 && //new
					!electron->isConv && //New
					electron->nMissingHits <= 0 && //new
					fabs(electron->dEtaIn) < 0.004 &&//elSCDeltaEtaEB &&
					fabs(electron->dPhiIn) < 0.03 && //elSCDeltaPhiEB &&
					electron->sieie < 0.01 && //elSigmaIEtaIEtaEB &&//SigmaIEtaIeta < elSigmaIEtaIEtaEB &&
          electron->hovere < 0.12 && //elHadOverEmEB  &&
					fabs(electron->d0) < 0.02 &&//elDxyEB &&
					fabs(electron->dz)  < 0.1 &&//elDzEB &&
					sumPfIso < elSumPFIsoEB 
				)      
      {                                                                                 
        elecList.push_back(electron);
				//std::cout << electron->pt << " " << electron->pfPt << std::endl;                                           
      } 
    }

		else if(fabs(electron->scEta) > 1.566){// > 1.479//electron->IsEE() == true){

      if(	fabs(electron->eta) < elEta && 
					electron->pt > elPt &&
					energyInversion < 0.05 && //new
					!electron->isConv && //New
					electron->nMissingHits <= 1 && //new
					fabs(electron->dEtaIn) < 0.007 && //elSCDeltaEtaEE &&//fabs(electron->SCDeltaeta) < elSCDeltaEtaEE &&
          fabs(electron->dPhiIn) < 0.03 &&//elSCDeltaPhiEE &&//fabs(electron->SCDeltaPhi()) < elSCDeltaPhiEE &&
					electron->sieie < 0.03 && //elSigmaIEtaIEtaEE &&                                                                 
          electron->hovere  < 0.10 &&// elHadOverEmEE  &&
        	fabs(electron->d0) < 0.02 &&//elDxyEE &&
        	fabs(electron->dz)  < 0.1 &&//elDzEE &&
					sumPfIso < elSumPFIsoEE )
      {
        elecList.push_back(electron);
      }
    }
    else ;//printf("Electron in gap\n");
	}	
}



void electronSelection(std::vector<baconhep::TGenParticle*> &elecList, TClonesArray* electrons, int id, std::vector<baconhep::TElectron*> &recoElecList ){
	if(id == 0) return;

	bool wasTau = false;
	bool push_back = false;

  for(Int_t i = 0; i < electrons->GetSize(); i++){
    baconhep::TGenParticle* electron = (baconhep::TGenParticle*) electrons->At(i);

		wasTau = false;
		if( electron == NULL ) break;
		if (abs(electron->pdgId) == 11 && electron->status == 1){ 

			if(electron->parent < 0 ) break;
      if( electrons->At( electron->parent ) == NULL) break;// return;
      baconhep::TGenParticle* mother;
      mother = (baconhep::TGenParticle*) electrons->At( electron->parent );
      baconhep::TGenParticle* a;

      while(abs(mother->pdgId) == 13 || abs(mother->pdgId) == 11 || abs(mother->pdgId) == 15){

				if( abs(mother->pdgId) == 15 ) wasTau = true;

				if(mother->parent < 0 ) break;
        if( electrons->At( mother->parent ) == NULL) break;//return;
        a = (baconhep::TGenParticle*)( electrons->At( mother->parent ) );
        mother = a;
      }

			if( abs(mother->pdgId) == id &&
					wasTau == true && 
					fabs(electron->eta) < elEta &&
					//fabs(electron->Dz((baconhep::TVertex*)vtxList.at(0)))  < muDz &&
					//fabs(electron->Dxy((baconhep::TVertex*)vtxList.at(0))) < muDxy &&
					electron->pt > elPt )
			{

				push_back = false;
				for( unsigned int j = 0; j < recoElecList.size(); j++ ){
					baconhep::TElectron* recoElectron = (baconhep::TElectron*) recoElecList.at(j);
					if( dR( electron->phi , recoElectron->phi, electron->eta, electron->eta) < 0.3 ){
						push_back = true;
					}
				}
				if( push_back == true){
					elecList.push_back((baconhep::TGenParticle*)electron);
					printf("electron %f", electron->pt);
				}
			}
		}
  }

}



