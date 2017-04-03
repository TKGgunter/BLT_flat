#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include <cstdlib> // for getenv

//FIXME: implement how to get trigger objects

using namespace baconhep;
TriggerSelector::TriggerSelector() {
	std::string cmssw_base = getenv("CMSSW_BASE");
	std::string trigfilename = cmssw_base + "/src/";
	trigfilename += "BaconAna/DataFormats/data/HLTFile_v2";


//	TTrigger some_ptr( "some_string");//some_string );
//	some_ptr.getTriggerBit("some string");
	_ttrigger.reset(new baconhep::TTrigger(trigfilename));
}

bool TriggerSelector::PassTrigger(const TriggerBits &iTrig) const {

	bool pass = true;
/*
	for(unsigned int iter = 0; iter < iTrig.size() ; iter++ ){
	
		if( iTrig[iter] == 1)	 printf( "iter %i\n", iter);
	}
*/
	for (const auto& triggerName: _triggerNames) {
		pass = _ttrigger->pass(triggerName, iTrig);
	}
	return pass;
}
