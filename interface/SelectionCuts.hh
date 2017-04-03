#ifndef SelectionCuts_hh
#define SelectionCuts_hh

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TVertex.hh"

#include "BLT/BLTAnalysis/interface/tgUtility.hh"

#include <TClonesArray.h>
#include <cmath>


//==============Good Muon Selection Cuts ==========
bool muIsGood									= true; 
float muEta										= 2.4;
float muNormChi2							= 10;
int muNumbValidHits						= 0;
float muSumPFIso							= 0.12;
float muPt										= 25.; // 20
bool muIsGLB									=	true;
bool muIsPF										=	true;
int muNumbMatchedStat				=	1;
int muTrackLayersMeasurements = 5;
float muDz										=	0.5;
float muDxy										=	0.2;


//==============Good Electron Selection Cuts=========
//bool elIsGood								= true;
//EB
float elSCDeltaEtaEB					= 0.004;
float elSCDeltaPhiEB					= 0.03;//0.06
float elSigmaIEtaIEtaEB				= 0.01;
float elHadOverEmEB						= 0.12;//0.01
float elSumPFIsoEB						= 0.1; //0.15
float elDzEB									= 0.1;
float elDxyEB									= 0.02;

//EE
float elSCDeltaEtaEE          = 0.0066; 
float elSCDeltaPhiEE          = 0.03;
float elSigmaIEtaIEtaEE       = 0.03;   
float elHadOverEmEE           = 0.10; //0.01; 
float elSumPFIsoEE            = 0.1; //0.15
float elDzEE									= 0.1;
float elDxyEE									= 0.02;

//General
float elPt										=	15.;
float elEta										= 2.5;

//==============Good Jet Selection Cuts==============
float jetNeuHadFrac						=	0.99;
float jetNeuEmFrac						= 0.99;
float jetNumConstit						=	1.;
float jetChEmFrac							=	.99;
float jetChHadFrac						=	0;
float jetPt										=	30.;	//Suggest cut is 25
float jetDR_mu_el							=	0.4;
float jetDR_me								= 0.3;
float jetEta									=	4.7;

float vtxCut									=	.3;

//void jetSelection( std::vector<TGenJet*> &jetList, TClonesArray* jets, std::vector<TMuon*> &muonList , std::vector<TElectron*> &elecList, TClonesArray* genParticles, int id);
void jetSelection( std::vector<baconhep::TJet*> &jetList , TClonesArray* jets, std::vector<baconhep::TMuon*> &muonList, std::vector<baconhep::TElectron*> &elecList, int nVtx);

bool test_bits_(unsigned int bits, unsigned int test);
void muonSelection( std::vector<baconhep::TMuon*> &muonList, TClonesArray* muons);
void muonSelection( std::vector<baconhep::TGenParticle*> &muonList, TClonesArray* muons, int x, std::vector<baconhep::TMuon*> &recoMuonList );
void softMuonSelection( std::vector<baconhep::TMuon*> &muonList, TClonesArray* muons);


void electronSelection( std::vector<baconhep::TElectron*> &elecList, TClonesArray* electrons, float rhoFactor);
void electronSelection( std::vector<baconhep::TGenParticle*> &elecList, TClonesArray* electrons, int x, std::vector<baconhep::TElectron*> &recoElecList);

#endif
