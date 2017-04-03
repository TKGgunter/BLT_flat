/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

// c++ libraries
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;

class WeightUtils: public TObject 
	{
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection, bool isRealData);

        void    Initialize() {};
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        float   GetPUWeight(float);
        //pair<float, float>   GetTriggerEffWeight(string, TLorentzVector&) const;
        float   GetMuonRecoEff(TLorentzVector&) const; 
        float   GetElectronRecoEff(TLorentzVector&) const; 

        ClassDef(WeightUtils, 0);

    private:
        //input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        TH1D*  puReweight;
        //TGraphErrors *_sf_IsoMu24_Eta2p1_data[3], *_sf_IsoMu24_Eta2p1_mc[3];

        TGraphErrors *_muSF2012_ID[4]; 
				TH2F* _elSF2012;
};

#endif

//#if !defined(__CINT__)
//CLassImp(WeightUtils);
//#endif
