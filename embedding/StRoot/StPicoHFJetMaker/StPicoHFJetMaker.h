#ifndef StPicoHFMyAnaMaker_h
#define StPicoHFMyAnaMaker_h

#include "StarClassLibrary/StThreeVectorF.hh"
#include "StarClassLibrary/StThreeVectorD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/SystemOfUnits.h"

#include "StPicoCuts/StPicoCuts.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoEmcTrigger.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StPicoEvent/StPicoMcTrack.h"
#include "StPicoEvent/StPicoMcVertex.h"

#include "../StPicoJetMaker/StPicoJetMaker.h"
#include "../StRefMultCorr/StRefMultCorr.h"

#include "TNtuple.h"
#include "TString.h"
#include "TSystem.h"

#include "TH2D.h"

#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include <vector>

#include "StMaker.h"

using namespace std;


class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StEmcADCtoEMaker;
class StBemcTables;


class StMcEvent;
class StMcTrack;
class StEvent;
class StTrack;
class StEmcDecoder;


class StPicoHFJetMaker : public StPicoJetMaker {
 public:
  //  StPicoHFJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName);
  //Jana changed to
  StPicoHFJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree);
  virtual ~StPicoHFJetMaker();
  
  virtual Int_t InitJets();
  virtual Int_t MakeJets();
  virtual void  ClearJets(Option_t *opt);
  virtual Int_t FinishJets();
	virtual Double_t GetTowerCalibEnergy(Int_t TowerId);
	virtual Double_t vertexCorrectedEta(double eta, double vz);
	virtual Bool_t GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx);
	virtual Int_t FindTriggerTowers(Int_t level);

    int SetHighTowerVar(StMcTrack* mcTrack, bool isele);


    StEmcADCtoEMaker *mADCtoEMaker;
  StBemcTables     *mTables;

  //void setRefMutCorr(StRefMultCorr *gRefMultCorr) { mRefmultCorrUtil = gRefMultCorr; }
  //StRefMultCorr* getRefMultCorr() { return mRefmultCorrUtil; }

    void doAuAu(bool kDoAuAu);
		void setEmbPythia(bool kEmbPythia);

    void setR(vector<float> &fR);
    void setAcuts(vector<float> &fAcuts);
    void setEmbPt(vector<float> &fEmbPt);

	void setGhostMaxrap(float fGhostMaxrap);
    void setNpTlead(int npTlead);
	void setR_bg(float fR_bg);
    void setNJetsRemove(int nJetsRemove);
    void setJetPtMin(float jetPtMin);
    void setCutETmin(float min);
    void setMcJetType(unsigned int us);
    
    void setTriggerThreshold(float fTrgthresh);
    
	/*void setRefMutCorr(StRefMultCorr* gRefMultCorr);
	StRefMultCorr* getRefMultCorr();*/
	void setRefMultCorr(StRefMultCorr* RefMultCorr);
	StRefMultCorr* getRefMultCorr();

    unsigned int mcJetType();
		
		
		void setHadronCorr(float corr);
		
		void setMaxNeutralFraction(float max);
		void setMaxDcaZHadronCorr(float max);
		
	void setMCparameters(float pThatmin, float pThatmax, float xweight);	

	enum eMcJetType {SingleParticle, Pythia, Jewel};

//    int InitRun(int runNumber);

 protected:
	TString   mInputFileName;        //! *.list - MuDst or picoDst  //Jana

private:

  // -- private members --------------------------


  // -- ADD USER MEMBERS HERE -------------------
    bool kDoAuAu;

	int fRunNumber;

	vector<float> fR;
	vector<float> fAcuts;
	vector<float> fEmbPt;

	bool kEmbPythia;

	float fRBg;
	float fGhostMaxrap;
	float fJetPtMin;
	int npTlead;
	int nJetsRemove;
	unsigned int mMcJetType;
	StRefMultCorr* mRefmultCorrUtil;
 	
 	// hadronic correction fraction
 	float fHadronCorr;
 	
	// -- max neutral fraction of a jet
	float maxneutralfrac;
	
	// -- max DCAz for global tracks used for hadronic correction 
	float	maxdcazhadroncorr;
	
	float fETmincut;
	
	const double mBarrelRadius = 225.405;
	
	//trigger threshold
	float fTrgthresh;
	//pThat range and cross-section weight
	float fpThatmin;
	float fpThatmax;
	float fWeight;

    float bemcEnergy[4801];
    int bemcADC[4801];
    StEmcDecoder* mEmcDecoder;
    StBemcTables* mBemcTables;
    StEvent* mEvent;


    // -- ADD USER MEMBERS HERE -------------------

  ClassDef(StPicoHFJetMaker, 0)

};


inline void StPicoHFJetMaker::doAuAu(bool kDoAuAu) {
	StPicoHFJetMaker::kDoAuAu = kDoAuAu;
}

inline void StPicoHFJetMaker::setEmbPythia(bool kEmbPythia){
	StPicoHFJetMaker::kEmbPythia = kEmbPythia;
}

inline void StPicoHFJetMaker::setR(vector<float> &fR) {
	StPicoHFJetMaker::fR = fR;
}

inline void StPicoHFJetMaker::setAcuts(vector<float> &fAcuts) {
	StPicoHFJetMaker::fAcuts = fAcuts;
}

inline void StPicoHFJetMaker::setEmbPt(vector<float> &fEmbPt) {
	StPicoHFJetMaker::fEmbPt = fEmbPt;
}

inline void StPicoHFJetMaker::setR_bg(float fR_bg) {
	StPicoHFJetMaker::fRBg = fR_bg;
}

inline void StPicoHFJetMaker::setGhostMaxrap(float fGhostMaxrap) {
	StPicoHFJetMaker::fGhostMaxrap = fGhostMaxrap;
}

inline void StPicoHFJetMaker::setNJetsRemove(int nJetsRemove) {
	StPicoHFJetMaker::nJetsRemove = nJetsRemove;
}

inline void StPicoHFJetMaker::setJetPtMin(float jetPtMin) {
	StPicoHFJetMaker::fJetPtMin = jetPtMin;
}

inline void StPicoHFJetMaker::setMcJetType(unsigned int us) {
	StPicoHFJetMaker::mMcJetType = us;
}

inline unsigned int StPicoHFJetMaker::mcJetType() {
	return mMcJetType;
}

inline void StPicoHFJetMaker::setNpTlead(int npTlead) {
    StPicoHFJetMaker::npTlead = npTlead;
}

inline void StPicoHFJetMaker::setHadronCorr(float corr)						{ fHadronCorr = corr;}

inline void StPicoHFJetMaker::setMaxNeutralFraction(float max)						{ maxneutralfrac = max;}

inline void StPicoHFJetMaker::setMaxDcaZHadronCorr(float max)						{ maxdcazhadroncorr = max;}

inline void StPicoHFJetMaker::setCutETmin(float min)						{ fETmincut = min;}
/*
inline void StPicoHFJetMaker::setRefMutCorr(StRefMultCorr *gRefMultCorr) {
	StPicoHFJetMaker::mRefmultCorrUtil = gRefMultCorr;
}*/

inline void StPicoHFJetMaker::setTriggerThreshold(float trgthresh) {fTrgthresh = trgthresh;}

inline void StPicoHFJetMaker::setRefMultCorr(StRefMultCorr *RefMultCorr) {
	StPicoHFJetMaker::mRefmultCorrUtil = RefMultCorr;
}

inline void StPicoHFJetMaker::setMCparameters(float pThatmin, float pThatmax, float xweight) {fpThatmin = pThatmin; fpThatmax = pThatmax; fWeight = xweight;}	

inline StRefMultCorr* StPicoHFJetMaker::getRefMultCorr() {
	return mRefmultCorrUtil;
}

#endif

//TODO: Add BEMC and TOF flag matching in the code.
//TODO: Add RefMultCorr routine to the current class.
