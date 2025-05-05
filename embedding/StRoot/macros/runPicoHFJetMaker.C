
/* **************************************************
 *   Run StPicoHFMyAnaMaker in different modes
 * --------------------------------------------------
 * run as :
 *  root -l -b -q StRoot/macros/loadSharedHFLibraries.C
 * StRoot/macros/runPicoHFMyAnaMaker.C++ or root -l -b -q
 * StRoot/macros/runPicoHFMyAnaMaker.C
 *
 * --------------------------------------------------
 *  - Different modes to use the  class
 *    - StPicoJetMaker::kAnalyze - don't write candidate trees, just fill
 * histograms inputFile : fileList of PicoDst files or single picoDst file
 *        outputFile: baseName for outfile
 *    - StPicoJetMaker::kWrite   - write candidate trees
 *        inputFile : path to single picoDist file
 *        outputFile: baseName for outfile
 *    - StPicoJetMaker::kRead    - read candidate trees and fill histograms
 *        inputFile : fileList of PicoDst files
 *        outputFile: baseName for outfile
 *
 * --------------------------------------------------
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)
 *
 * **************************************************
 */

#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

#include "StChain/StChain.h"
#include "StChain/StMaker.h"

#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "St_db_Maker/St_db_Maker.h"

#include "StPicoCuts/StPicoCuts.h"

#include "loadSharedHFLibraries.C"

#include <cstdio>
#include <ctime>
#include <iostream>

#include "StPicoHFJetMaker/StPicoHFJetMaker.h"

#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

using namespace std;

#else
class StChain;
#endif

StChain *chain;

void runPicoHFJetMaker(
    TString inputFile, TString outputFile = "outputBaseName",
    const unsigned int makerMode = 0 /*kAnalyze*/, TString treeName = "picoDst",
    const float pThatmin = 0., float pThatmax = -1.,
    const float xweight = 1.) { // this line added for Embedding analysis

#ifdef __CINT__
  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  chain = new StChain();
  // ========================================================================================
  // makerMode    = StPicoJetMaker::kAnalyze;
  // ========================================================================================

  cout << "Maker Mode    " << makerMode << endl;
  cout << "TreeName      " << treeName << endl;
  cout << "Input file: " << inputFile << endl;
  cout << "Loading bad run list: " << endl;
  TString badRunListFileName = "BadRunList_14.list";
  // check if exists
  if (!gSystem->AccessPathName(badRunListFileName)) {
    cout << "Bad run list file found!" << endl;
  } else {
    cout << "Bad run list file not found! Exiting..." << endl;
    exit(1);
  }

  if (makerMode == StPicoJetMaker::kAnalyze) {
    if (!inputFile.Contains(".list") && !inputFile.Contains("picoDst.root")) {
      cout << "No input list or picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  } else if (makerMode == StPicoJetMaker::kWrite) {
    if (!inputFile.Contains("picoDst.root")) {
      cout << "No input picoDst root file provided! Exiting..." << endl;
      exit(1);
    }
  } else {
    cout << "Unknown makerMode! Exiting..." << endl;
    exit(1);
  }

  // cout <<  "inputFile.Data() is  " << inputFile.Data() << endl;

  StMessMgr *msg = StMessMgr::Instance();
  msg->SwitchOff("Could not make BEMC detector");

  StPicoDstMaker *picoDstMaker =
      new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile,
                         "picoDstMaker"); // for local testing only
  St_db_Maker *dbMaker = new St_db_Maker("StarDb", "MySQL:StarDb");
  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
  StPicoHFJetMaker *stPicoHFJetMaker =
      new StPicoHFJetMaker("stPicoHFJetMaker", picoDstMaker, outputFile, "");

  stPicoHFJetMaker->setMakerMode(makerMode);
  stPicoHFJetMaker->setTreeName(treeName);
  // stPicoHFJetMaker->setMcMode(true);
  stPicoHFJetMaker->setMcMode(false);

  StPicoCuts *picoCuts = new StPicoCuts("PicoCuts");
  stPicoHFJetMaker->setPicoCuts(picoCuts);

  // ---------------------------------------------------
  // -- Set Base cuts for HF analysis

  // -- File name of bad run list
  picoCuts->setBadRunListFileName(badRunListFileName);

  // -- File name of hot tower list
  // picoCuts->setHotTowerListFileName(hotTowerListFileName);

  // -- ADD USER CUTS HERE ----------------------------
  // picoCuts->setCutVzMax(6.); //HFT range
  picoCuts->setCutVzMax(30.);
  picoCuts->setCutVzVpdVzMax(3.);
  // picoCuts->setCutRefMult(396, 100000); //REMEMBER now only central
  picoCuts->setCutRefMult(0, 100000);

  // 2014 MB HFT triggers

  /*picoCuts->addTriggerId(450050);    // vpdmb-5-p-nobsmd-hlt
  picoCuts->addTriggerId(450060);    // vpdmb-5-p-nobsmd-hlt
  picoCuts->addTriggerId(450005);    // vpdmb-5-p-nobsmd
  picoCuts->addTriggerId(450015);    // vpdmb-5-p-nobsmd
  picoCuts->addTriggerId(450025);    // vpdmb-5-p-nobsmd*/

  // 2014 MB triggers
  picoCuts->addTriggerId(450010); // vpdmb-30
  picoCuts->addTriggerId(450020); // vpdmb-30
  picoCuts->addTriggerId(450008); // vpdmb-5
  picoCuts->addTriggerId(450018); // vpdmb-5

  // 2014 BHT2*VPD30 triggers
  //	picoCuts->addTriggerId(450202);    // BHT2*VPD30
  //	picoCuts->addTriggerId(450212);    // BHT2*VPD30

  // 2014 BHT3 triggers
  // picoCuts->addTriggerId(450203);    // BHT3
  // picoCuts->addTriggerId(450213);    // BHT3

  // starprod AuAu 2014 luminosity - split low mid high, without lumi = early
  // runs

  /*
          //SL16j triggers
    picoCuts->addTriggerId(520802);    // VPDMB-5-p-hlt
    picoCuts->addTriggerId(520812);    // VPDMB-5-p-hlt
    picoCuts->addTriggerId(520822);    // VPDMB-5-p-hlt
    picoCuts->addTriggerId(520832);    // VPDMB-5-p-hlt
    picoCuts->addTriggerId(520842);    // VPDMB-5-p-hlt

    picoCuts->addTriggerId(520001);    // VPDMB-5-p-sst
    picoCuts->addTriggerId(520011);    // VPDMB-5-p-sst
    picoCuts->addTriggerId(520021);    // VPDMB-5-p-sst
    picoCuts->addTriggerId(520031);    // VPDMB-5-p-sst
    picoCuts->addTriggerId(520041);    // VPDMB-5-p-sst
    picoCuts->addTriggerId(520051);    // VPDMB-5-p-sst
  */
  picoCuts->setCutNHitsFitMin(14);
  picoCuts->setCutNHitsFitnHitsMax(0.52);
  picoCuts->setCutRequireHFT(false);
  picoCuts->setCutTPCNSigma(30.0);

  picoCuts->setCutDcaMin(1.0);
  picoCuts->setCutEta(1);
  picoCuts->setCutPtRange(0.2, 30.0); // default
  picoCuts->setCutERange(0.2, 30.0);  //

  vector<float> R;
  R.push_back(0.2);
  R.push_back(0.3);
  R.push_back(0.4);

  vector<float> Acuts;
  Acuts.push_back(0.07);
  Acuts.push_back(0.2);
  Acuts.push_back(0.4);

  vector<float> EmbPt;
  EmbPt.push_back(3.0);
  EmbPt.push_back(5.0);
  EmbPt.push_back(7.0);
  EmbPt.push_back(10.0);
  EmbPt.push_back(20.0);

  // TPC setters




  stPicoHFJetMaker->setNpTlead(10);
  stPicoHFJetMaker->setGhostMaxrap(1.0);
  stPicoHFJetMaker->setR(R);
  stPicoHFJetMaker->setJetPtMin(0.2); // default
  // stPicoHFJetMaker->setJetPtMin(0.5);
  stPicoHFJetMaker->setCutETmin(0.2);
  stPicoHFJetMaker->setAcuts(Acuts);
  stPicoHFJetMaker->setNJetsRemove(1);
  stPicoHFJetMaker->setR_bg(0.3);
  stPicoHFJetMaker->setEmbPt(EmbPt);

  stPicoHFJetMaker->setEmbPythia(0); // 0 = single particle, 1 = pythia6
  // stPicoHFJetMaker->setEmbPythia(1);

  stPicoHFJetMaker->setHadronCorr(1.0);
  // stPicoHFJetMaker->setHadronCorr(0.);

  stPicoHFJetMaker->setTriggerThreshold(
      18); //~4.2 GeV is the HT2 threshold = 18 ADC

  stPicoHFJetMaker->setMaxNeutralFraction(0.95); // default
  // stPicoHFJetMaker->setMaxNeutralFraction(95); //turn off for neutral jets

  stPicoHFJetMaker->setMaxDcaZHadronCorr(
      3.0); // cm, max DCA_z for global tracks used for hadronic correction

  stPicoHFJetMaker->setMCparameters(pThatmin, pThatmax,
                                    xweight); // pThat range and xsection weight
  // set refmultCorr
  // cout<<"test"<<endl;

  // USE gRefMultCorr or RefMultCorr??

  // StRefMultCorr* grefmultCorrUtil =
  // CentralityMaker::instance()->getgRefMultCorr_P18ih(); //new StRefMultCorr,
  // info about Run16, SL16d in the same file as for Run14, SL16d

  StRefMultCorr *grefmultCorrUtil =
      CentralityMaker::instance()
          ->getRefMultCorr_P18ih_VpdMB30_MidLow(); // new StRefMultCorr, info
                                                   // about Run14
  stPicoHFJetMaker->setRefMultCorr(grefmultCorrUtil);

  // StRefMultCorr* refmultCorrUtil =
  // CentralityMaker::instance()->getgRefMultCorr_VpdMB30(); //new
  // StRefMultCorr, info about Run14
  // stPicoHFJetMaker->setRefMultCorr(refmultCorrUtil);
  // cout<<"test2"<<endl;
  //  ========================================================================================

  // ========================================================================================

  // start_t = clock(); // getting starting time
  chain->Init();
  // cout << "chain->Init();" << endl;
  long int nEvents = picoDstMaker->chain()->GetEntries();
  //  int nEvents = 10000;
  // cout << " Total entries = " << nEvents << endl;
  // if(nEvents>total) nEvents = total;
  for (Int_t i = 0; i < nEvents; i++) {
    // if(i%10000 == 0) cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) {
      cout << "Bad return code!" << iret << endl;
      break;
    }
    // cout << "Evt. no. " << i << endl;
    // total++;
  }
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!" << endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  // double duration = (double) (clock() - start_t) / (double)CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;
  // cout << "Time needed " << duration << " s" << endl;
  cout << "****************************************** " << endl;

  delete chain;
}
