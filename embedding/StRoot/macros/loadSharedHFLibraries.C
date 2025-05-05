#include "TSystem.h"

void loadSharedHFLibraries() {

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  TString fastjetPath = gSystem->Getenv("FASTJET");
  if (!fastjetPath.IsNull()) {
      gSystem->Load(Form("%s/lib/libfastjet", fastjetPath.Data()));
  } else {
      std::cerr << "FASTJET environment variable is not set!" << std::endl;
  }
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical"); 
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

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
	//PYTHIA6 libraries
	gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");
	//gSystem->Load("StarPythia6.so");	

  cout << " loading of shared HF libraries are done" << endl;

  // -->>> ADD your own library/class HERE 
	gSystem->Load("StPicoHFJetMaker"); //analysis maker should be last
 }
