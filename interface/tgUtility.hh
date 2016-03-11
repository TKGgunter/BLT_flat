#ifndef tgUtility_h
#define tgUtility_h

#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TPhoton.hh"
#include "BaconAna/DataFormats/interface/TJet.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BLT/BLTAnalysis/interface/tgPhysObject.hh"
#include <TVector2.h>
#include <iostream>
#include <string>
#include <list>


void testFunc(std::string testStr, float testNumb = -99 );

//=======================================
float qtCalc(baconhep::TMuon* lep1, baconhep::TMuon* lep2 );

float qtCalc(baconhep::TMuon* lep1, baconhep::TElectron* lep2 );

float qtCalc(baconhep::TElectron* lep1, baconhep::TElectron* lep2 );

float qtCalc(baconhep::TElectron* lep1, baconhep::TMuon* lep2 );

float qtCalc(TGPhysObject* lep1, TGPhysObject* lep2);
//=========================================

float deltaPhi(float dPhi_);

//float deltaPhi(baconhep::TMuon* lep1, baconhep::TElectron* lep2);			//*******

//float deltaPhi(baconhep::TMuon* lep1, baconhep::TElectron* lep2, TJet* jet);		//******

float dR(baconhep::TMuon* lep1, baconhep::TJet* jet);

float dR(baconhep::TElectron* lep1, baconhep::TJet* jet);   
/*
float dR(baconhep::TMuon* lep1, TGenJet* jet);

float dR(baconhep::TElectron* lep1, TGenJet* jet);
*/
/////float dR(TGenParticle* gParticle, TGenJet* jet);

/////float dR(TJet* jet, TGenJet* genJet);
/*
float dR(baconhep::TMuon* lep1, TMET* met);

float dR(baconhep::TElectron* lep1, TMET* met);

float dotLepMet(baconhep::TMuon* lep, TMET* met );

float dotLepMet(baconhep::TElectron* lep, TMET* met);
*/
//=====================================

void tgSort(std::vector<baconhep::TMuon*>& muonList);

void tgSort(std::vector<baconhep::TElectron*>& elecList);

void tgSort(std::vector<baconhep::TJet*>& jetList);

//////void tgSort(std::vector<TGenJet*>& jetList);

void tgSort(std::vector<TGPhysObject*>& objList);

//=======================================

void tgCleanVector(std::vector<baconhep::TMuon*>& muonList);

void tgCleanVector(std::vector<baconhep::TElectron*>& elecList);

void tgCleanVector(std::vector<baconhep::TJet*>& jetList, int nVtx);

void tgConcateList(std::vector<baconhep::TMuon*> &muonList, std::vector<baconhep::TElectron*> &elecList, std::vector<TGPhysObject*> &objList);
#endif
