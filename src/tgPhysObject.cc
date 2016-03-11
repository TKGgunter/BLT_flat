#include "BLT/BLTAnalysis/interface/tgPhysObject.hh"

/*
    TGPhysObject();
    TGPhysObject(TMuon* obj);
    TGPhysObject(TElectron* obj);
    TGPhysObject(TJet* obj);

    ~TGPhysObject(){}

    float pt, eta, phi;
    int id;

    float Mag();
*/
TGPhysObject::TGPhysObject(){
	pt	= 0.;
	eta	= 0.;
	phi = 0.;

	p4.SetPtEtaPhiM(pt, eta, phi, 0);
  id	= -999;
  ELE_MASS = 0.000511;
  MUON_MASS = 0.105658369;
	

}

TGPhysObject::TGPhysObject(baconhep::TMuon* muon){
  pt  = muon->pfPt; 
  eta = muon->pfEta; 
  phi = muon->pfPhi;

	p4.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  id  = 13;
}

TGPhysObject::TGPhysObject(baconhep::TElectron* electron){
  pt  = electron->pfPt; 
  eta = electron->pfEta; 
  phi = electron->pfPhi;

	p4.SetPtEtaPhiM(pt, eta, phi, ELE_MASS);
  id  = 11;
}

TGPhysObject::TGPhysObject(baconhep::TJet* jet){
  pt  = jet->pt;
  eta = jet->eta;
  phi = jet->phi;

	p4.SetPtEtaPhiM(pt, eta, phi, 1);
  id  = jet->mcFlavor;
}
//float TGPhysObject::Mag(){
//	
//
//}

