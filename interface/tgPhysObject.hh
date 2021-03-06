#ifndef TGPHYSOBJECTS_HH
#define TGPHYSOBJECTS_HH

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include <TObject.h>
#include <TLorentzVector.h>

class TGPhysObject{
	public:
	TGPhysObject();	
	TGPhysObject(baconhep::TMuon* muon);
	TGPhysObject(baconhep::TElectron* electron);
	TGPhysObject(baconhep::TJet* jet);
	TGPhysObject(baconhep::TGenParticle* particle);

	~TGPhysObject(){}

	void calcP4();
	void decontructP4();

	double ELE_MASS; // = 0.000511;
	double MUON_MASS;// = 0.105658369;

	float pt, eta, phi, iso;
	int id, charge, mother;
	TLorentzVector p4;

};
#endif  
