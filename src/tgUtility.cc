#include "BLT/BLTAnalysis/interface/tgUtility.hh"
#include <cmath>

void testFunc(std::string testStr, float testNumb){
	if( testNumb == -99) std::cout << "Is " << testStr << " working." << std::endl;

	else std::cout << "Is " << testStr << " working. Are you expecting the folling value: " << testNumb << std::endl; 
}

//======================================

float qtCalc(baconhep::TMuon* lep1, baconhep::TMuon* lep2 ){
	float qt = 0;
	qt = sqrt(pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2) + pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2));
	return qt; 
}

float qtCalc(baconhep::TMuon* lep1, baconhep::TElectron* lep2 ){
	float qt = 0;
	qt = sqrt(pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2) + pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2)); 
	return qt;
}

float qtCalc(baconhep::TElectron* lep1, baconhep::TElectron* lep2 ){
	float qt = 0;
	qt = sqrt(pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2) + pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2)); 
	return qt;
}

float qtCalc(baconhep::TElectron* lep1, baconhep::TMuon* lep2 ){
	float qt = 0;
	qt = sqrt(pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2) + pow(lep1->pfPt*sin(lep1->pfPhi) + lep2->pfPt*sin(lep2->pfPhi), 2)); 
	return qt;
}

float qtCalc(TGPhysObject* lep1, TGPhysObject* lep2){
	float qt = 0;
	qt = sqrt(pow(lep1->pt*sin(lep1->phi) + lep2->pt*sin(lep2->phi),2) + pow(lep1->pt*cos(lep1->phi) + lep2->pt*cos(lep2->phi), 2));
	return qt;

}

//=====================================

float deltaPhi(float dPhi){
	if(dPhi != dPhi){
		printf("Phi%f\n", dPhi);
		dPhi = 0;
	}
	dPhi = TVector2::Phi_mpi_pi(dPhi);
	//if( dPhi > M_PI) dPhi = abs(dPhi - M_PI);
	
	return dPhi;
}

float dR(baconhep::TMuon* lep1, baconhep::TJet* jet){
	float dR_ = 0.;
	float dPhi = deltaPhi(lep1->pfPhi-jet->phi);
	float dEta = lep1->pfEta - jet->eta;
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR(baconhep::TElectron* lep1, baconhep::TJet* jet){  
	float dR_ = 0.; 
	float dPhi = deltaPhi(lep1->pfPhi-jet->phi);   
	float dEta = lep1->pfEta - jet->eta; 
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR(TGPhysObject* obj1, TGPhysObject* obj2){  
	float dR_ = 0.; 
	float dPhi = deltaPhi(obj1->phi - obj2->phi);   
	float dEta = obj1->eta - obj2->eta; 
	dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
	return dR_;
}

float dR( float phi1, float phi2, float eta1, float eta2){
	float dR_  = 0.;
	float dPhi = deltaPhi( phi1-phi2);
	float dEta = eta1-eta2;
	dR_ = sqrt( pow(dPhi, 2.) + pow(dEta, 2) ); 
	return dR_;
}

/*
float dR(baconhep::TMuon* lep1, TGenJet* jet){                                                     
  float dR_ = 0.;
  float dPhi = deltaPhi(lep1->pfPhi-jet->pfPhi);                                        
  float dEta = lep1->pfEta - jet->pfEta;
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                              
  return dR_;                                                                           
}                                                                                       

float dR(baconhep::TElectron* lep1, TGenJet* jet){                                                 
  float dR_ = 0.; 
  float dPhi = deltaPhi(lep1->pfPhi-jet->pfPhi);                                        
  float dEta = lep1->pfEta - jet->pfEta; 
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                              
  return dR_;                                                                           
}

float dR(TGenParticle* gParticle, TGenJet* jet){
  float dR_ = 0.;
  float dPhi = deltaPhi(gParticle->pfPhi-jet->pfPhi);
  float dEta = gParticle->pfEta - jet->pfEta;
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));
  return dR_;
}

float dR(TJet* jet, TGenJet* genJet){
  float dR_ = 0.;
  float dPhi = deltaPhi(genJet->pfPhi-jet->pfPhi);                                    
  float dEta = genJet->pfEta - jet->pfEta;                                            
  dR_ = sqrt(pow(dPhi, 2) + pow(dEta, 2));                                               
  return dR_;  
}
*/
/*
float dR(baconhep::TMuon* lep1, TMET* met){                                                     
  float dR_ = 0.;
	//dR_ = met->DeltaR(*lep1);
  return dR_;                                                                           
}                                                                                       

float dR(baconhep::TElectron* lep1, TMET* met){                                                 
  float dR_ = 0.; 
  //dR_ = met->DeltaR(*lep1);                                        
  return dR_;                                                                           
}                                                                                       
*/
//===================================
/*
float dotLepMet(baconhep::TMuon* lep, TMET* met ){
	float dotProduct = 0;
	dotProduct = lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi) + lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi);
	return dotProduct;
}

float dotLepMet(baconhep::TElectron* lep, TMET* met){
	float dotProduct = 0;
	dotProduct = lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi) + lep->pfPt*sin(lep->pfPhi) * met->pfPt*sin(met->pfPhi);  
	return dotProduct;
}
*/  
//====================================


void tgSort(std::vector<baconhep::TMuon*> & muonList){
	float pfPt_1 = 0.;
	float pfPt_2 = 0.;
	for(unsigned int i=0; i < muonList.size(); i++){

		baconhep::TMuon* muon_1 = (baconhep::TMuon*) muonList.at(i);
		if (i == 0) continue;
		
		for(unsigned int j = 0; j < i ; j++ ){
			baconhep::TMuon* muon_2 = (baconhep::TMuon*) muonList.at(j);  
			pfPt_1 = muon_1->pfPt;
			pfPt_2 = muon_2->pfPt;
	//		printf("%f \t %f\n", pfPt_2, pfPt_1);
			if(pfPt_2 < pfPt_1 ){
				muonList.insert(muonList.begin()+j, muon_1);
				muonList.erase(muonList.begin()+i+1);	
				continue;
			}
		}
	}
/*
  for(unsigned int i= 0; i < muonList.size(); i++){
    baconhep::TMuon* muon = (baconhep::TMuon*) muonList.at(i);
    printf("muonPT:\t%f\t", muon->pfPt);
  }
  printf("\n");
*/
}

void tgSort(std::vector<baconhep::TElectron*> & elecList){
  float pfPt_1 = 0.;
  float pfPt_2 = 0.;
  for(unsigned int i=0; i < elecList.size(); i++){
    
    baconhep::TElectron* elec_1 = (baconhep::TElectron*) elecList.at(i);                                          
    if (i == 0) continue;                                                             
    for(unsigned int j = 0; j < i ; j++ ){
      baconhep::TElectron* elec_2 = (baconhep::TElectron*) elecList.at(j);                                        
      pfPt_1 = elec_1->pfPt;                                                              
      pfPt_2 = elec_2->pfPt;                                                              
      //printf("%f \t %f\n", pfPt_2, pfPt_1);                                                 
      if(pfPt_2 < pfPt_1 ){                                                                 
        elecList.insert(elecList.begin()+j, elec_1);                                    
        elecList.erase(elecList.begin()+i+1); 
        continue;
      }
    }                                                                                   
  } 
}

void tgSort(std::vector<baconhep::TJet*> & jetList){
  float pfPt_1 = 0.;
  float pfPt_2 = 0.;
  for(unsigned int i=0; i < jetList.size(); i++){
    
    baconhep::TJet* jet_1 = (baconhep::TJet*) jetList.at(i);                                          
    if (i == 0) continue;                                     
    for(unsigned int j = 0; j < i ; j++ ){
      baconhep::TJet* jet_2 = (baconhep::TJet*) jetList.at(j);                                        
      pfPt_1 = jet_1->pt;                                                              
      pfPt_2 = jet_2->pt;                                                              
      if(pfPt_2 < pfPt_1 ){                                                                 
        jetList.insert(jetList.begin()+j, jet_1);                                    
        jetList.erase(jetList.begin()+i+1); 
        continue;
      }
    }                                                                                   
  } 
}


void tgSort(std::vector<TGPhysObject*> & objList){
  float pfPt_1 = 0.;
  float pfPt_2 = 0.;
  for(unsigned int i=0; i < objList.size(); i++){

    TGPhysObject* obj_1 = (TGPhysObject*) objList.at(i);
    if (i == 0) continue;
    for(unsigned int j = 0; j < i ; j++ ){
      TGPhysObject* obj_2 = (TGPhysObject*) objList.at(j);
      pfPt_1 = obj_1->pt;
      pfPt_2 = obj_2->pt;
      if(pfPt_2 < pfPt_1 ){
        objList.insert(objList.begin()+j, obj_1);
        objList.erase(objList.begin()+i+1);
        continue;
      }
    }
  }
}


/*
void tgSort(std::vector<TGenJet*>& jetList){
  float pfPt_1 = 0.;                                                                       
  float pfPt_2 = 0.; 
  for(unsigned int i=0; i < jetList.size(); i++){                                        
    
    TGenJet* jet_1 = (TGenJet*) jetList.at(i);                                               
    if (i == 0) continue;                                                                
    for(unsigned int j = 0; j < i ; j++ ){
      TGenJet* jet_2 = (TGenJet*) jetList.at(j);                                             
      pfPt_1 = jet_1->pfPt;                                                                
      pfPt_2 = jet_2->pfPt;                                                                
      if(pfPt_2 < pfPt_1 ){                                                                  
        jetList.insert(jetList.begin()+j, jet_1);                                        
        jetList.erase(jetList.begin()+i+1);                                              
        continue;                                                                        
      }                                                                                  
    }                                                                                    
  } 	
}
*/
//================================================
void tgCleanVector(std::vector<baconhep::TMuon*>& muonList){
	for(unsigned int i = 1 ; i < muonList.size(); i++ ){
		baconhep::TMuon* muon1 = (baconhep::TMuon*) muonList.at(i);
		baconhep::TMuon* muon2 = (baconhep::TMuon*) muonList.at(i-1);
		if (muon1->pfPt == muon2->pfPt && muon1->pfPhi == muon2->pfPt){ 
			muonList.erase(muonList.begin() + i);
			i=i-1;
		}
	}
}

void tgCleanVector(std::vector<baconhep::TElectron*>& elecList){
  for(unsigned int i = 1 ; i < elecList.size(); i++ ){
    baconhep::TElectron* electron1 = (baconhep::TElectron*) elecList.at(i);
    baconhep::TElectron* electron2 = (baconhep::TElectron*) elecList.at(i-1);
    if (electron1->pfPt == electron2->pfPt && electron1->pfPhi == electron2->pfPhi){
      elecList.erase(elecList.begin() + i);
			i=i-1;
    }
  }
}


void tgCleanVector(std::vector<baconhep::TJet*>& jetList){
  for(unsigned int i = 1 ; i < jetList.size(); i++ ){
    baconhep::TJet* jet1 = (baconhep::TJet*) jetList.at(i);
    baconhep::TJet* jet2 = (baconhep::TJet*) jetList.at(i-1);
		if( jet1 == NULL || jet2 == NULL)
		{
			printf("jets broken");
			break;
		}
    if (jet1->pt == jet2->pt){
      jetList.erase(jetList.begin() + i);
      i=i-1;
    }
		//printf("nVtx %i\n", nVtx);
		/*
		if(fabs(jet2->eta) < 2.5){
			if ( jet2->betaStar/log(nVtx-0.64) > .2  && jet2->dR2Mean > 0.06) {
				printf( "kill 1 %i\n", i );
				jetList.erase( jetList.begin() + i - 1);
			}
		}
		else if(fabs(jet2->eta) < 2.75  ){
			if (jet2->betaStar/log(nVtx-0.64) > .3  && jet2->dR2Mean > 0.05){ 
				printf( "kill 2\n" );
				jetList.erase( jetList.begin() + i -1  );
			}
		}
		else if( fabs( jet2->eta) < 3){
			if ( jet2->dR2Mean > 0.05) { 
				printf( "kill 3\n" );
				jetList.erase( jetList.begin() + i -1  ); 
			}
		}

		else {
			if( jet2->dR2Mean > 0.055){
				printf("kill 4\n");
				jetList.erase(jetList.begin() + i - 1 );
			}
		}
		*/
  }
}


//=================================================
//void tgDeltaR( std::vector<baconhep::TElectron*> electList, std::vector<baconhep::TMuon*> muonList){
//
//}

//=================================================

void tgConcateList(std::vector<baconhep::TMuon*> &muonList, std::vector<baconhep::TElectron*> &elecList, std::vector<TGPhysObject*>& objList){

	for(unsigned int i = 0; i < elecList.size() ; i++){
		baconhep::TElectron* electron = (baconhep::TElectron*) elecList.at(i);
		TGPhysObject* po_ele = new TGPhysObject(electron);
		objList.push_back(po_ele);
	}
	for(unsigned int i = 0; i < muonList.size() ; i++){
		baconhep::TMuon* muon = (baconhep::TMuon*) muonList.at(i);
		TGPhysObject* po_muon = new TGPhysObject(muon);
		objList.push_back(po_muon);
	}

	tgSort(objList);	
}

void tgConcateList(std::vector<baconhep::TGenParticle*> &muonList, std::vector<baconhep::TGenParticle*> &elecList, std::vector<TGPhysObject*>& objList){

	for(unsigned int i = 0; i < elecList.size() ; i++){
		baconhep::TGenParticle* electron = (baconhep::TGenParticle*) elecList.at(i);
		TGPhysObject* po_ele = new TGPhysObject(electron);
		objList.push_back(po_ele);
	}
	for(unsigned int i = 0; i < muonList.size() ; i++){
		baconhep::TGenParticle* muon = (baconhep::TGenParticle*) muonList.at(i);
		TGPhysObject* po_muon = new TGPhysObject(muon);
		objList.push_back(po_muon);
	}

	tgSort(objList);	
}

