#include "BLT/BLTAnalysis/interface/SelectionCuts.hh"
//#include "tgUtility.C"

bool test_bits_(unsigned int bits, unsigned int test) {
    return (bits & test) == test;
}
/*
void jetSelection( std::vector<TGenJet*> &jetList , TClonesArray* jets, std::vector<TMuon*> &muonList, std::vector<TElectron*> &elecList, TClonesArray* genParticles , int id ){

  bool isJetGood = true;
  for(Int_t i = 0; i < jets->GetSize(); i++){
    TGenJet* jet = (TGenJet*) jets->At(i);
    isJetGood = true;
    if( (jet->pfEta < jetEta) &&
    //    (jet->NeuHadFrac() < jetNeuHadFrac) &&
    //    (jet->NeuEmFrac() < jetNeuEmFrac) &&
    //    (jet->NumConstit() > jetNumConstit ) &&  //Commented out for genJets
    //    (jet->ChEmFrac() < jetChEmFrac) &&
    //    (jet->ChHadFrac() > jetChHadFrac) &&
        jet->pfPt > jetPt )//&&                                                            
      //  jet->NumChPart() > 0 )									//Commented  out for genJets             
    {                                                                                   
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
        int count = 0;     //1 for gen particles                                                          
        for( int j = 0; j < genParticles->GetSize() ; j++){  
            TGenParticle* gJet = (TGenParticle*) genParticles->At(j);
				//		printf("test\n");
						float deltaR = dR(gJet, jet);
				//		if ( deltaR < .4) std::cout << deltaR << " " << abs(gJet->GetPDGId()) << " "  << gJet->GetStatus() <<" pt genPart: " << gJet->pfPt << " pT Jet: " << jet->pfPt  <<std::endl;
						if ( deltaR < .4 && (abs(gJet->GetPDGId()) < 11 || (abs(gJet->GetPDGId()) > 21 && abs(gJet->GetPDGId()) < 24)) ){//&& gJet->GetStatus() == 1){
//							printf("test\n");
							if(gJet->Mother() != NULL){
								TGenParticle* mother;
								mother  = gJet->Mother();
								TGenParticle a;

								while(abs(mother->GetPDGId()) < 6 || abs(mother->GetPDGId()) == 21 || abs(mother->GetPDGId()) == 22  ){

									if(mother->Mother() == NULL) break;//return;
									a = *(mother->Mother());
									mother = &a;
								}

								if (abs(mother->GetPDGId()) == 24 || abs(mother->GetPDGId()) == 23 || abs(mother->GetPDGId()) == 5 || abs(mother->GetPDGId()) == 6) count++;
						}
						else if(abs(gJet->GetPDGId()) == 24 || abs(gJet->GetPDGId()) == 23 || abs(gJet->GetPDGId()) == 22 || abs(gJet->GetPDGId()) == 6 || abs(gJet->GetPDGId()) == 5 ) count++;
					}
				} 
		//		printf("/=====/\n");
				if(count > 0){
					//printf("asdf %i \n",jet->JetFlavor() );
					jetList.push_back( (TGenJet*)jet);
				}
			}
		}
	}
}
*/

//Muon Selection
void muonSelection(std::vector<baconhep::TMuon*> &muonList, TClonesArray* muons){//, TClonesArray* primaryVtx ){
/*  std::vector<TPrimaryVtx*> vtxList;

  for(int i = 0; i < primaryVtx->GetSize() ; i++){
    TPrimaryVtx* primVtx = (TPrimaryVtx*) primaryVtx->At(i);
    if( !primVtx->IsFake() && primVtx->NDof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ){
			vtxList.push_back(primVtx);
		}
	}
*/
	for(Int_t i = 0; i < muons->GetEntries(); i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muons->At(i);
		float sumPfIso = (muon->chHadIso03 + muon->neuHadIso03 + muon->gammaIso03) / muon->pfPt;//(muon->PfIsoChargedHad() + muon->PfIsoNeutral()+muon->PfIsoPhoton())/muon->pfPt;

		if (
				fabs(muon->pfEta) < muEta   &&
				//muon->IsGood() == muIsGood &&
				test_bits_(muon->typeBits, baconhep::kGlobal) == muIsGLB && //muon->IsGLB() == muIsGLB &&
				test_bits_(muon->typeBits, baconhep::kPFMuon) == muIsPF &&//muon->IsPF() == muIsPF	&&
				int(muon->nMatchStn) >  muNumbMatchedStat &&
				int(muon->nTkLayers)  > muTrackLayersMeasurements &&//muon->TrackLayersWithMeasurement() > muTrackLayersMeasurements &&
				muon->muNchi2 < muNormChi2 &&
				muon->nValidHits > muNumbValidHits &&
				fabs(muon->dz) < muDz &&//Dz((TPrimaryVtx*)vtxList.at(0)))  < muDz &&
				fabs(muon->d0) < muDxy &&//Dxy((TPrimaryVtx*)vtxList.at(0))) < muDxy &&
				sumPfIso < muSumPFIso && 
				muon->pfPt > muPt)
		{
			
			muonList.push_back(muon);
		}
	}
} 



/*
void muonSelection(std::vector<TMuon*> &muonList, TClonesArray* muons, int id, TClonesArray* primaryVtx ){

  std::vector<TPrimaryVtx*> vtxList;

  for(int i = 0; i < primaryVtx->GetSize() ; i++){
    TPrimaryVtx* primVtx = (TPrimaryVtx*) primaryVtx->At(i);
    if( !primVtx->IsFake() && primVtx->NDof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ) vtxList.push_back(primVtx);
  }


	if (id == 0) return;
	 
  for(Int_t i = 0; i < muons->GetSize(); i++){
    TGenParticle* muon = (TGenParticle*) muons->At(i);


		if ( abs(muon->GetPDGId()) == 13 && muon->GetStatus() == 1){
			if(muon->Mother() == NULL) break;//return;
			TGenParticle* mother;
			mother  = muon->Mother();
			TGenParticle a;

			while(abs(mother->GetPDGId()) == 13 || abs(mother->GetPDGId()) == 11 || abs(mother->GetPDGId()) == 15){

				if(mother->Mother() == NULL) break;//return;
				a = *(mother->Mother());
				mother = &a;
			}
			
			//if (abs(mother->GetPDGId()) != id) break;// != id return; 

			if ( abs(mother->GetPDGId()) == id &&
					fabs(muon->pfEta) < muEta   &&
		      fabs(muon->Dz((TPrimaryVtx*)vtxList.at(0)))  < muDz &&
	        fabs(muon->Dxy((TPrimaryVtx*)vtxList.at(0))) < muDxy &&
					muon->pfPt > muPt)
			{
				muonList.push_back((TMuon*)muon);
			}
		}
  }
}
*/




//Electron Selection
void electronSelection(std::vector<baconhep::TElectron*> &elecList, TClonesArray* electrons){//, TClonesArray* primaryVtx){
/*
  std::vector<TPrimaryVtx*> vtxList;

  for( int i = 0; i < primaryVtx->GetSize() ; i++){
    TPrimaryVtx* primVtx = (TPrimaryVtx*) primaryVtx->At(i);
    if( !primVtx->IsFake() && primVtx->NDof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ) vtxList.push_back(primVtx);
  }
*/
	for(Int_t i = 0; i < electrons->GetEntries(); i++){
		baconhep::TElectron* electron = (baconhep::TElectron*) electrons->At(i);
		float sumPfIso = (electron->chHadIso03 + electron->neuHadIso03 + electron->gammaIso03 )/electron->pfPt;
		
		if(fabs(electron->pfEta) < 1.479 ){//electron->IsEB() == true){

			if(	fabs(electron->pfEta) < elEta &&
					electron->pfPt > elPt &&
					fabs(electron->dEtaIn) < elSCDeltaEtaEB &&
					fabs(electron->dPhiIn) < elSCDeltaPhiEB &&
					electron->sieie < elSigmaIEtaIEtaEB &&//SigmaIEtaIeta < elSigmaIEtaIEtaEB &&
          electron->hovere < elHadOverEmEB  &&
					fabs(electron->dz)  < elDzEB &&
					fabs(electron->d0) < elDxyEB )//&&
					 //sumPfIso < elSumPFIsoEB )      
      {                                                                                 
        elecList.push_back(electron);                                                   
      } 
    }

		else if(fabs(electron->pfEta) > 1.479){//electron->IsEE() == true){

      if(	fabs(electron->pfEta) < elEta && 
					electron->pfPt > elPt &&
					fabs(electron->dEtaIn) < elSCDeltaEtaEE &&//fabs(electron->SCDeltaeta) < elSCDeltaEtaEE &&
          fabs(electron->dPhiIn) < elSCDeltaPhiEE &&//fabs(electron->SCDeltaPhi()) < elSCDeltaPhiEE &&
					electron->sieie < elSigmaIEtaIEtaEE &&                                                                 
          electron->hovere  < elHadOverEmEE  &&
        	fabs(electron->dz)  < elDzEE &&
        	fabs(electron->d0) < elDxyEE &&
					sumPfIso < elSumPFIsoEE )
      {
        elecList.push_back(electron);
      }
    }
    else printf("Electron in gap\n");
	}	
}


/*
void electronSelection(std::vector<TElectron*> &elecList, TClonesArray* electrons, int id, TClonesArray* primaryVtx ){
	if(id == 0) return;

	std::vector<TPrimaryVtx*> vtxList;

  for( int i = 0; i < primaryVtx->GetSize() ; i++){
    TPrimaryVtx* primVtx = (TPrimaryVtx*) primaryVtx->At(i);
    if( !primVtx->IsFake() && primVtx->NDof() > 4. && fabs(primVtx->z()) <= 24. && fabs(primVtx->Perp()) <=2  ) vtxList.push_back(primVtx);
  }

  for(Int_t i = 0; i < electrons->GetSize(); i++){
    TGenParticle* electron = (TGenParticle*) electrons->At(i);

		if (abs(electron->GetPDGId()) == 11 && electron->GetStatus() == 1){ 

      if(electron->Mother() == NULL) break;// return;
      TGenParticle* mother;
      mother  = electron->Mother();
      TGenParticle a;

      while(abs(mother->GetPDGId()) == 13 || abs(mother->GetPDGId()) == 11 || abs(mother->GetPDGId()) == 15){

        if(mother->Mother() == NULL) break;//return;
        a = *(mother->Mother());
        mother = &a;
      }

			//printf("id: %i\n", mother->GetPDGId());
//      if (abs(mother->GetPDGId()) != id) break;//return;

			if( abs(mother->GetPDGId()) == id &&
					fabs(electron->pfEta) < elEta &&
					fabs(electron->Dz((TPrimaryVtx*)vtxList.at(0)))  < muDz &&
					fabs(electron->Dxy((TPrimaryVtx*)vtxList.at(0))) < muDxy &&
					electron->pfPt > elPt )
			{
				elecList.push_back((TElectron*)electron);
			}
		}
  }
}


*/
