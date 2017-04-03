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
  mother = 0;
	iso = 0;
}

TGPhysObject::TGPhysObject(baconhep::TMuon* muon){
  pt  = muon->pfPt; 
  eta = muon->pfEta; 
  phi = muon->pfPhi;
	charge = muon->q;

	p4.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  id  = 13;
	mother = 0;
	iso = muon->trkIso03;
}

TGPhysObject::TGPhysObject(baconhep::TElectron* electron){
  pt  = electron->pt;//pfPt; 
  eta = electron->eta;//pfEta; 
  phi = electron->phi;//pfPhi;
	charge = electron->q;

	p4.SetPtEtaPhiM(pt, eta, phi, ELE_MASS);
  id  = 11;
	mother = 0;
	iso = electron->chHadIso03 + electron->neuHadIso03 + electron->gammaIso03;
}

TGPhysObject::TGPhysObject(baconhep::TJet* jet){
  pt  = jet->pt;
  eta = jet->eta;
  phi = jet->phi;
	charge = jet->q;

	p4.SetPtEtaPhiM(pt, eta, phi, 1);
  id  = jet->mcFlavor;
	mother = 0;
	iso = 0;
}

TGPhysObject::TGPhysObject(baconhep::TGenParticle* particle){
  pt  = particle->pt; 
  eta = particle->eta; 
  phi = particle->phi;
	charge = 0;

	p4.SetPtEtaPhiM(pt, eta, phi, MUON_MASS);
  id  = -999;
	mother = 0;
	iso = 0;
}

void TGPhysObject::calcP4(){
	p4.SetPtEtaPhiM(pt, eta, phi, 0);
}

void TGPhysObject::decontructP4(){
	pt = p4.Pt();
	eta = p4.Eta();
	phi = p4.Phi();

}



