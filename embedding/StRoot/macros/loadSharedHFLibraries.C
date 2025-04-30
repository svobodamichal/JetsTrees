#include "TSystem.h"

#ifndef __CINT__
#include "StMuDSTMaker/COMMON/macros/loadSharedLibraries.C"
#endif

extern TSystem* gSystem;

void loadSharedHFLibraries() {

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("fastjet1/fastjet_install/lib/libfastjet");
  gSystem->Load("fastjet1/fastjet_install/lib/libsiscone");
  gSystem->Load("fastjet1/fastjet_install/lib/libsiscone_spherical"); 
  gSystem->Load("fastjet1/fastjet_install/lib/libfastjetplugins");
  gSystem->Load("fastjet1/fastjet_install/lib/libfastjettools");
  gSystem->Load("fastjet1/fastjet_install/lib/libfastjetcontribfragile");
 	gSystem->Load("StBTofUtil");
	gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoCuts");
  gSystem->Load("StPicoJetMaker");
  gSystem->Load("StRefMultCorr"); 
  gSystem->Load("St_db_Maker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); //rember, order of loading makers matters
	gSystem->Load("StPicoHFJetMaker"); //analysis maker should be last

	//PYTHIA6 libraries
	gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");	

  cout << " loading of shared HF libraries are done" << endl;

  // -->>> ADD your own library/class HERE 

 }
