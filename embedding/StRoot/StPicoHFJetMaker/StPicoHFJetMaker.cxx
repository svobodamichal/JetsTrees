#include "StPicoHFJetMaker.h"
#include "JetInfo.h"

#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"

#include "BemcNewCalib.h"

//FastJet 3
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
/*
#include <../../fastjet/config.h>
#include <../../fastjet/PseudoJet.hh>
#include <../../fastjet/JetDefinition.hh>
#include <../../fastjet/ClusterSequence.hh>
#include <../../fastjet/ClusterSequenceArea.hh>
//#ifdef FASTJET_VERSION
#include <../../fastjet/Selector.hh>
#include <../../fastjet/tools/Subtractor.hh>
#include <../../fastjet/tools/JetMedianBackgroundEstimator.hh>*/

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
//#ifdef FASTJET_VERSION
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

//#include </gpfs01/star/pwg/licenrob/jets/fastjet/contrib/SoftDrop.hh>	//robotmon 
//#include </gpfs01/star/pwg/licenrob/jets/fastjet/contrib/RecursiveSoftDrop.hh>
//#include <fastjet/internal/base.hh> //robotmon
//#endif
#include <vector>

#include <TRandom3.h>

using namespace std;  //robotmon
using namespace fastjet;
//using namespace contrib;	//robotmon

ClassImp(StPicoHFJetMaker)

bool trackErr = false;
bool towErrPlus = false;
bool towErrMinus = false;



//_________________________________match MC tracks to TPC tracks...has to be here, fastjet cannot be included in StPicoHFJetMaker.h file___________________________________________________________
//not used
bool MatchTracks(vector<PseudoJet> &McTracks,vector<PseudoJet> &RcTracks){
	bool found = false;
	for (unsigned int i = 0; i < McTracks.size(); i++){
		double mceta = McTracks[i].eta();	
		double mcphi = McTracks[i].phi();	
		double mcpt = McTracks[i].pt();	
		double cut = 0; //delta r cut
		if (mcpt < 1) cut = 0.0475; //cut from deltar distribution minimum - pT dependent
		else if (mcpt < 3) cut = 0.0285;
		else cut = 0.014;
		int nmatch = 0; //number of matched tracks
		int ind = 0; //index of matched track
		vector<int> jvec; //indices of matched track candidates

		for (unsigned int j = 0; j < RcTracks.size(); j++){
			if (RcTracks[j].user_index() > 1000000) continue; //skip already-matched tracks
			double rceta = RcTracks[j].eta();
			double rcphi = RcTracks[j].phi();	
			double deltaeta = mceta-rceta;
			double deltaphi = mcphi-rcphi;
			if (deltaphi > TMath::Pi()) deltaphi = 2.0*TMath::Pi() - deltaphi;
			if (deltaphi < -1.0*TMath::Pi()) deltaphi = 2.0*TMath::Pi() + deltaphi;
			double deltar = sqrt(pow(deltaeta,2)+pow(deltaphi,2));
			if (deltar < cut){ //cout << mcpt << " track match with deltar " << deltar << endl;
					nmatch++;
					jvec.push_back(j);		
					found = true;
			}
		}
		if (nmatch == 0) continue;
		if (nmatch > 1) {//cout << nmatch <<  " matches" << endl; //if more than 1 match, select random match
				srand (time(NULL)); //random seed
				ind = rand() % nmatch; //integer between 0 and nmatch-1 (hopefully)
				//cout << nmatch << " random match number " << ind << endl;
			}

		McTracks[i].set_user_index(1000000+i);
		//cout << McTracks[i].user_index() << endl;
		RcTracks[jvec[ind]].set_user_index(1000000+i);
		//cout << "next mc track" << endl;
	}


	return found;
}

//_________________________________rewritten jet matching using matched tracks using idTruth and qaTruth...has to be here, fastjet cannot be included in StEffEmbed.h file___________________________________________________________
bool MatchJets(vector<PseudoJet> McJets, vector<PseudoJet> Rcjets, vector<double> McPtLeads, vector<double> Rcleads, vector<pair<PseudoJet,PseudoJet>>* matched, vector<pair<double,double>>* matchedPtLead,
		vector<pair<double,double>>* matchedNeutralFraction, double R){
	bool found = false;
	vector<PseudoJet> RcJets = Rcjets; //copy RC jets, so we can remove after match
	vector<double> RcPtLeads = Rcleads; //copy RC leading particles, so we can remove after match
	vector<double> matchpT,matchR,matchpTfrac,matchR2,matchpTfrac2,matchdeltapT; //must be cleared at the end (of first loop)
	vector<pair<PseudoJet,PseudoJet>> matchedTmp, matchedTmp2;	//must be cleared at the end
	vector<pair<double,double>> matchedPtLeadTmp;	//must be cleared at the end
	vector<pair<double,double>> matchedNeutralFractionTmp;//must be cleared at the end
	vector<int> jvec; //indices of matched jet candidates
	vector<int> mcindex; //user indices of MC tracks, must be cleared at the end
	vector<double> matchtrackpT; //pT from matched tracks, must be cleared at the end

    for (unsigned int i = 0; i < McJets.size(); i++) {
		found = false;
		vector<PseudoJet> constituentsMc = sorted_by_pt(McJets[i].constituents());
		double nfractionMc = 0;
		double neutralpTMc = 0;
		double pT_jetMc = McJets[i].perp();
//		if(pT_jetMc>10.2 && pT_jetMc<10.22) {
//            cout << "MC pT  " << i << "  " << pT_jetMc << "  " << McJets[i].eta() << "  " << McJets[i].phi() << endl;

            for (unsigned int ic = 0; ic < constituentsMc.size(); ++ic){
//                if (constituentsMc[ic].perp() > 0.1) {cout<<"MC const " << ic << "    " << constituentsMc[ic].perp() <<"  " <<constituentsMc[ic].eta()<<"   "<<constituentsMc[ic].phi()<<" "<< constituentsMc[ic].user_index()  <<endl;}

                int uidxMC = constituentsMc[ic].user_index();
                if (uidxMC > -1) mcindex.push_back(uidxMC);//select matched mc tracks
                if (uidxMC == 0) neutralpTMc += constituentsMc[ic].perp();
            }
            nfractionMc = neutralpTMc / pT_jetMc;
//            cout << "Počet recontructed  " << RcJets.size() << endl;

            for (unsigned int j = 0; j < RcJets.size(); j++) {
                double pT_jetRc = RcJets[j].perp();
                vector <PseudoJet> constituentsRc = sorted_by_pt(RcJets[j].constituents());
                double nfractionRc = 0;
                double neutralpTRc = 0;
                double pTmatch = 0;
//                cout << "RC pT  " << j << "  " << pT_jetRc << "  " << RcJets[j].eta() << "  " << RcJets[j].phi()<< endl;


                double etaMC = McJets[i].eta();
                double phiMC = McJets[i].phi();
                double etaRC = RcJets[j].eta();
                double phiRC = RcJets[j].phi();
                double etaDiff = etaMC - etaRC;
                double phiDiff = phiMC - phiRC;
                if(phiDiff<=-TMath::Pi()){phiDiff=phiDiff+TMath::TwoPi();}
                if(phiDiff>=TMath::Pi()){phiDiff=phiDiff-TMath::TwoPi();}


                for (unsigned int irc = 0; irc < constituentsRc.size(); ++irc) {
                    int uidx = constituentsRc[irc].user_index();
                    double constpT = constituentsRc[irc].perp();
//                    if (constpT > 0.1) { cout << irc << "    " << constpT <<"  " <<constituentsRc[irc].eta()<<"   "<<constituentsRc[irc].phi()<<"  "<<constituentsRc[irc].user_index()  <<endl; }

                    std::vector<int>::iterator it;
                    it = std::find(mcindex.begin(), mcindex.end(), uidx);
                    int idx = std::distance(mcindex.begin(), it);
                    if (uidx > 0 && it != mcindex.end()) {
                        pTmatch += constpT; //search for RC track user index in the mcidx vector
                    } else if (uidx == 0 || uidx == 9999) { neutralpTRc += constpT; }//select towers
                }
                if (pTmatch > 0) {
                    found = true; // found at least one RC track which matches track at MC level
                    matchtrackpT.push_back(pTmatch);
                    nfractionRc = neutralpTRc / pT_jetRc;
                    matchedNeutralFractionTmp.push_back(make_pair(nfractionMc, nfractionRc));
                    matchedTmp.push_back(make_pair(McJets[i], RcJets[j]));
                    matchedPtLeadTmp.push_back(make_pair(McPtLeads[i], RcPtLeads[j]));


                    jvec.push_back(j);
                } //end of if (pTmatch > 0) - track match candidate found
            } //end of RC jets loop
            if (!found) {
                jvec.clear();
                mcindex.clear();
                matchtrackpT.clear();
                matchedTmp.clear();
                matchedPtLeadTmp.clear();
                matchedNeutralFractionTmp.clear();
                continue;
            } //no match

            std::vector<double>::iterator result;
            result = std::max_element(matchtrackpT.begin(), matchtrackpT.end());
            int index = std::distance(matchtrackpT.begin(), result);

            //apply area and eta cut
            RcJets.erase(RcJets.begin() + jvec[index]); //remove already-matched det-lvl jet
            RcPtLeads.erase(RcPtLeads.begin() + jvec[index]); //and corresponding pTlead
            //cout << "removed jet no. " << jvec[index] << endl;
            double area = matchedTmp[index].second.area();
            double Area_cuts[3] = {0.07, 0.2, 0.4}; //stupid way to "access" fAcuts
            int ac = R * 10 - 2; //find index R=0.2 -> 0, R=0.3 -> 1, R=0.4 -> 2

            if (area > Area_cuts[ac] && fabs(matchedTmp[index].second.eta()) < 1 - R &&
                fabs(matchedTmp[index].first.eta()) < 1 - R && matchedNeutralFractionTmp[index].second < 0.95) {
                matched->push_back((pair <PseudoJet, PseudoJet>) matchedTmp[index]);
                matchedPtLead->push_back(matchedPtLeadTmp[index]);
                matchedNeutralFraction->push_back(matchedNeutralFractionTmp[index]);
            }
            //}

            jvec.clear();
            mcindex.clear();
            matchtrackpT.clear();
            matchedTmp.clear();
            matchedPtLeadTmp.clear();
            matchedNeutralFractionTmp.clear();
//        } //test event

    	}

// for (unsigned int i = 0; i < matched->size(); i++) cout << matched->at(i).first.perp() << " " << matched->at(i).second.perp() << endl;


    return found;
}
bool MatchJetsEtaPhi(vector<PseudoJet> McJets, vector<PseudoJet> Rcjets, vector<double> McPtLeads, vector<double> Rcleads, vector<pair<PseudoJet, PseudoJet>>* matched, vector<pair<double, double>>* matchedPtLead, vector<pair<double, double>>* matchedNeutralFraction, double R, vector<double>* differEta, vector<double>* differPhi, vector<pair<double, double>>* deltaR) {
    bool found = false;
    vector<PseudoJet> RcJets = Rcjets; //copy RC jets, so we can remove after match
    vector<double> RcPtLeads = Rcleads; //copy RC leading particles, so we can remove after match
    vector<pair<PseudoJet, PseudoJet>> matchedTmp;    //must be cleared at the end
    vector<pair<double, double>> matchedPtLeadTmp;    //must be cleared at the end
    vector<pair<double, double>> matchedNeutralFractionTmp; //must be cleared at the end
    vector<int> jvec; //indices of matched jet candidates
    //vector<int> mcindex; //user indices of MC tracks, must be cleared at the end
    // vector<double> matchtrackpT; //pT from matched tracks, must be cleared at the end
    vector<double> rdiffer; //difference in r = sqrt(deltaeta^2+deltaphi^2)

    //cout << "MC jet pT, R = " << R << endl;


//    cout << "MC jets size " << McJets.size() << " RC jets size " << RcJets.size() << endl;
    for (unsigned int i = 0; i < McJets.size(); i++) {
        found = false;
        vector<PseudoJet> constituentsMc = sorted_by_pt(McJets[i].constituents());
        double nfractionMc = 0;
        double neutralpTMc = 0;
        const PseudoJet& mcJet = McJets[i];
        double mcEta = mcJet.eta();
        double mcPhi = mcJet.phi();
        double pT_jetMc = mcJet.perp();

        // cout << pT_jetMc << endl;

        for (unsigned int ic = 0; ic < constituentsMc.size(); ++ic){
            int uidxMC = constituentsMc[ic].user_index();
            //  if (uidxMC > -1) mcindex.push_back(uidxMC);//select matched mc tracks
            if (uidxMC == 0) neutralpTMc += constituentsMc[ic].perp();

        }
        nfractionMc = neutralpTMc / pT_jetMc;

     //   cout<<"MC pT "<< pT_jetMc<<" MC eta  " << mcEta << " MC phi  " << mcPhi <<endl;
        for (unsigned int j = 0; j < RcJets.size(); j++) {
            const PseudoJet& rcJet = RcJets[j];
            double rcEta = rcJet.eta();
            double rcPhi = rcJet.phi();
            double etaDiff = mcEta - rcEta;
            double phiDiff = mcPhi - rcPhi;
            if (phiDiff <= -TMath::Pi()) { phiDiff = phiDiff + TMath::TwoPi(); }
            if (phiDiff >= TMath::Pi()) { phiDiff = phiDiff - TMath::TwoPi(); }
            double deltar = sqrt(etaDiff*etaDiff+phiDiff*phiDiff);
            double pT_jetRc = rcJet.perp();
        //     cout<<"RC eta  " << rcEta << "RC phi  " << rcPhi <<" RC pT "<< pT_jetRc<<endl;

            //cout << "RC pT: " << pT_jetRc << " deltar " << deltar << endl;

            vector<PseudoJet> constituentsRc = sorted_by_pt(rcJet.constituents());
            double nfractionRc = 0;
            double neutralpTRc = 0;
            // double pTmatch = 0;
            differEta->push_back(etaDiff);
            differPhi->push_back(phiDiff);
            deltaR->push_back(make_pair(deltar, pT_jetRc));

            for (unsigned int irc = 0; irc < constituentsRc.size(); ++irc) {
                int uidx = constituentsRc[irc].user_index();
                double constpT = constituentsRc[irc].perp();

                //std::vector<int>::iterator it;
                // it = std::find(mcindex.begin(), mcindex.end(), uidx);
                // int idx = std::distance(mcindex.begin(), it);
                // if (uidx > 0 && it != mcindex.end()) {
                //     pTmatch += constpT; //search for RC track user index in the mcidx vector
                /*} else */if (uidx == 0 || uidx == 9999) { neutralpTRc += constpT; } //select towers
            }

            // Check if RC jet and MC jet match based on spatial properties
            //if (abs(etaDiff) < 0.1 && abs(phiDiff) < 0.12) { //= 0.156 when summed in quadrature
            if (deltar < 0.6*R) { //motivated by ALICE analysis, to be revised
                found = true;
                //matchtrackpT.push_back(pTmatch);
                rdiffer.push_back(deltar);
                nfractionRc = neutralpTRc / pT_jetRc;
                matchedNeutralFractionTmp.push_back(make_pair(nfractionMc, nfractionRc));
                matchedTmp.push_back(make_pair(mcJet, rcJet));
                matchedPtLeadTmp.push_back(make_pair(McPtLeads[i], RcPtLeads[j]));
                jvec.push_back(j);
            }
        } // RC loop

        if (!found) {
            jvec.clear();
            rdiffer.clear();
            //   mcindex.clear();
            // matchtrackpT.clear();
            matchedTmp.clear();
            matchedPtLeadTmp.clear();
            matchedNeutralFractionTmp.clear();
            continue;
        } //no match

        std::vector<double>::iterator result;
        // result = std::max_element(matchtrackpT.begin(), matchtrackpT.end());
        // int index = std::distance(matchtrackpT.begin(), result);
        result = std::min_element(rdiffer.begin(),rdiffer.end());
        int index = std::distance(rdiffer.begin(), result);
        //cout << "No. matches " << matchedTmp.size() << " index " << index << endl;
        //apply area and eta cut
        RcJets.erase(RcJets.begin() + jvec[index]); //remove already-matched det-lvl jet
        RcPtLeads.erase(RcPtLeads.begin() + jvec[index]); //and corresponding pTlead

        double area = matchedTmp[index].second.area();
        double Area_cuts[3] = { 0.07, 0.2, 0.4 };
        int ac = R * 10 - 2;

        //cout << "Matched pT MC " << matchedTmp[index].first.pt() << " RC "  << matchedTmp[index].second.pt() << endl;


        if (area > Area_cuts[ac] && fabs(matchedTmp[index].second.eta()) < 1 - R &&
            fabs(matchedTmp[index].first.eta()) < 1 - R && matchedNeutralFractionTmp[index].second < 0.95) {

            //cout << "PAIR ACCEPTED" << endl;

            matched->push_back(make_pair(matchedTmp[index].first, matchedTmp[index].second));
            matchedPtLead->push_back(make_pair(matchedPtLeadTmp[index].first, matchedPtLeadTmp[index].second));
            matchedNeutralFraction->push_back(make_pair(matchedNeutralFractionTmp[index].first, matchedNeutralFractionTmp[index].second));
        }

        /*else { if (area < Area_cuts[ac]) cout << "area cut failed" << endl; else if (fabs(matchedTmp[index].second.eta()) > 1-R) cout << "RC eta cut failed" << endl;else if (fabs(matchedTmp[index].first.eta()) > 1-R) cout << "MC eta cut failed" << endl; else if (matchedNeutralFractionTmp[index].second > 0.95) cout << "NEF cut failed" << endl;}*/

        jvec.clear();
        rdiffer.clear();
        //   mcindex.clear();
        //   matchtrackpT.clear();
        matchedTmp.clear();
        matchedPtLeadTmp.clear();
        matchedNeutralFractionTmp.clear();
    } // MC loop

    return found;
}

StPicoHFJetMaker::StPicoHFJetMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName, char const* inputHFListHFtree = "") :
StPicoJetMaker(name, picoMaker, outputBaseFileName, inputHFListHFtree),  mRefmultCorrUtil(NULL) {


  // constructor
}

// _________________________________________________________
StPicoHFJetMaker::~StPicoHFJetMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFJetMaker::InitJets() {


  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread")); 
	assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());

    int npTleadbins = 25;
    float pTleadmin = 0;
    float pTleadmax = 25;
    int nptbins = 1000;
    float ptminbin = -40;
    float ptmaxbin = 60;
    int netabins = 100*2;
    float etaminbin = -1;
    float etamaxbin = 1;
    int nphibins = 120;
    float phiminbin = 0;//-TMath::Pi();
    float phimaxbin = TMath::TwoPi();
    int zbins = 50;
    float zmin = -50;
    float zmax = 50;
    int npttrackbins = 120;
    float pttrackmin = 0;
    float pttrackmax = 30;
		
    int ncetabins = 100*2;
    float cetaminbin = -1;
    float cetamaxbin = 1;
    int ncphibins = 120;
    float cphiminbin = 0;//-TMath::Pi();
    float cphimaxbin = TMath::TwoPi();

		int nptembbins = 200;
    float ptembminbin = 0;
    float ptembmaxbin = 25;
    float deltaptembminbin = -30;
    float deltaptembmaxbin = 50;
    
    TH1::SetDefaultSumw2();

    mOutList->Add(new TH1D("hweight", "weight", 135, 0, 3));
	mOutList->Add(new TH1D("hcent", "centrality", 10, -1, 9));

    //General track QA
	mOutList->Add(new TH2D("heta_phi_tr", "track phi vs. eta;#phi [-];#eta [-]", nphibins, phiminbin, phimaxbin, netabins, etaminbin, etamaxbin));
	mOutList->Add(new TH1D("hphi_tr", "track phi;#phi", nphibins, phiminbin, phimaxbin));
	mOutList->Add(new TH1D("heta_tr", "track eta;#eta", netabins, etaminbin, etamaxbin));    
	mOutList->Add(new TH1D("hpT_tr", "track pT; p_{T} (GeV/c)", npttrackbins, pttrackmin, pttrackmax));
//	mOutList->Add(new TH1D("hE_tr", "track E from BEMC; E [GeV]", nEtrackbins, Etrackmin, Etrackmax));
//	mOutList->Add(new TH1D("hpTBEMC_tr", "track pT from BEMC; p_{T} (GeV/c)", nptBEMCtrackbins, ptBEMCtrackmin, ptBEMCtrackmax));
	mOutList->Add(new TH2D("hdca_z_tr", "track DCA vs z-vertex", 90, 0, 3, zbins , zmin, zmax));
    	mOutList->Add(new TH1D("hdca_tr", "track DCA", 90, 0, 3));
    	mOutList->Add(new TH2D("hdca_pT", "track DCA vs. p_{T}", 90, 0, 3, npttrackbins, pttrackmin, pttrackmax));
	mOutList->Add(new TH1D("hcharged_tr", "track charge", 90, 0, 3));	
	
	//MC tracks
	mOutList->Add(new TH2D("hMcEtaPhi", "MC particle eta vs phi; #eta (-); #phi (-)", netabins, etaminbin, etamaxbin,nphibins, phiminbin, phimaxbin));	
		
		

		//mOutList->Add(new TH1D("hdca_XY_tr", "global track DCA_XY; DCA_{xy} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH1D("hdca_Z_tr", "global track DCA_Z; DCA_{z} [cm]", 200, -20, 20));	
    //mOutList->Add(new TH2D("hdca_XYZ_tr", "global track DCA_Z vs DCA_XY; DCA_{z} [cm]; DCA_{xy} [cm]", 200, -20, 20, 200, -20, 20));	
   	//for (int lumi = 0; lumi < 5; lumi++) 
		//mOutList->Add(new TH2D(Form("hdca_XYZ_tr_lumi%i", lumi), Form("global track DCA_Z vs DCA_XY, %i < lumi < %i kHz; DCA_{z} [cm]; DCA_{xy} [cm]", lumi*15, (lumi+1)*15), 200, -20, 20, 200, -20, 20));	


		//mOutList->Add(new TH1D("hEcorr", "E from BTow after correction; E [GeV]", 340, -4.5, 29.5));

  //  mOutList->Add(new TH2D("hE_tow", "energy vs tower ID;tower ID;E [GeV]", 4800, 1, 4800, 300, -0.5, 29.5));
    //mOutList->Add(new TH2D("hp_tow", "track momentum vs tower ID;tower ID;p (GeV/c)", 4800, 1, 4800, 300, -0.5, 29.5));
  //mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));		
	  mOutList->Add(new TH2D("heta_phi_tow", "tower eta vs phi; #eta [-]; #phi [-]", netabins/5, etaminbin, etamaxbin, nphibins, phiminbin, phimaxbin));	
	  mOutList->Add(new TH1D("hET_tow", "tower ET; E_{T} (GeV)", npttrackbins, pttrackmin, pttrackmax));
      mOutList->Add(new TH1D("hADC", "tower ADC; Nevím jendotky", 100, 0, 100));


    if (isMcMode()) { //not used
    }

	TString hname ="hrho";
	//mOutList->Add(new TH1D(hname,"median energy density; #rho_{ch} [GeV/sr]",50,0,50));		
        hname = "hrho_mult";
        //mOutList->Add(new TH2D(hname,"median energy density vs multiplicity; RefMult [-] ; #rho_{ch} [GeV/sr]",600, 0, 600, 50,0,50));
        hname = "hfrho";
        mOutList->Add(new TH1D(hname,"median energy density; #rho [GeV/sr]",100,0,100));
        hname = "hfrho_mult";
        mOutList->Add(new TH2D(hname,"median energy density vs multiplicity ; RefMult [-] ; #rho [GeV/sr]",600, 0, 600, 100,0,100));

 	if (!isMcMode()) {
        for(unsigned int r = 0; r < fR.size(); r++) {

            hname = Form("hphi_MCRC_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "phi MC-RC", nphibins, phiminbin, phimaxbin));
            hname = Form("heta_MCRC_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "eta MC-RC", netabins, etaminbin, etamaxbin));
            hname = Form("hEtaPhi_MC-RC_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname, "MC-RC #eta, #phi; #eta (-); #phi (-)", netabins, etaminbin, etamaxbin,nphibins, phiminbin, phimaxbin));
            hname = Form("hphi_MCRCw_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "phi MC-RC", nphibins, phiminbin, phimaxbin));
            hname = Form("heta_MCRCw_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH1D(hname, "eta MC-RC", netabins, etaminbin, etamaxbin));
            hname = Form("hEtaPhi_MC-RCw_R0%.0lf",fR[r]*10);
            mOutList->Add(new TH2D(hname, "MC-RC #eta, #phi; #eta (-); #phi (-)", netabins, etaminbin, etamaxbin,nphibins, phiminbin, phimaxbin));
//            hname = Form("hDeltaR_R0%.0lf",fR[r]*10);
//            mOutList->Add(new TH2D(hname, "deltaR vs reco pT",nptbins, ptminbin, ptmaxbin,100, 0, 1));
//            hname = Form("hDeltaRw_R0%.0lf",fR[r]*10);
//            mOutList->Add(new TH2D(hname, "deltaR vs reco pT",nptbins, ptminbin, ptmaxbin,100, 0, 1));


            //TString hname = Form("hpT_pTlead_R0%.0lf",fR[r]*10);
            //CHARGED JETS
            //mOutList->Add(new TH2D(hname, "jet pTcorr vs pTleading; p_{T} (GeV/c); p_{T}^{lead} (GeV/c)", nptbins, ptminbin, ptmaxbin, npTleadbins, pTleadmin, pTleadmax));
	/*	hname = Form("hpT_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "pT; p_{T} (GeV/c)", nptbins, ptminbin, ptmaxbin)); 
		hname = Form("heta_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH1D(hname, "jet eta;#eta", netabins, etaminbin, etamaxbin));
	    	hname = Form("hphi_R0%.0lf",fR[r]*10);
            	mOutList->Add(new TH1D(hname, "jet phi;#phi", nphibins, phiminbin, phimaxbin));
          	hname = Form("hjetarea_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname,"jet area", 100, 0, 1));
	        hname = Form("hjetpTarea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"jet pTmeasured vs area; A [-]; p_{T} (GeV/c)" ,100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hjetpTcorrArea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"jet pTreco vs area; A [-]; p_{T}^{reco} (GeV/c)",100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hnparticlesinjet_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"#particles in jet vs jet pT; # of particles; p_{T}^{lead} (GeV/c)", 300*fR[r], 0, 300*fR[r], npTleadbins, pTleadmin, pTleadmax));
		hname = Form("hceta_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));
		hname = Form("hcphi_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH1D(hname, "jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));
*/
		//FULL JETS
		hname = Form("hfceta_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet constituents eta;#eta", ncetabins, cetaminbin, cetamaxbin));
		hname = Form("hfcphi_R0%.0lf",fR[r]*10);
            	mOutList->Add(new TH1D(hname, "full jet constituents phi;#phi", ncphibins, cphiminbin, cphimaxbin));	
		hname = Form("hfjetarea_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname,"full jet area", 100, 0, 1));
	        hname = Form("hfjetpTarea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"full jet pTmeasured vs area; A [-]; p_{T} (GeV/c)", 100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hfjetpTcorrArea_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"full jet pTreco vs area; A [-]; p_{T}^{reco} (GeV/c)", 100, 0, 1, nptbins, ptminbin, ptmaxbin));
	        hname = Form("hfnparticlesinjet_R0%.0lf",fR[r]*10);
        	mOutList->Add(new TH2D(hname,"#particles in full jet vs full jet pT; # of particles; p_{T}^{lead} (GeV/c)", 500*fR[r], 0, 500*fR[r], npTleadbins, pTleadmin, pTleadmax));
		hname = Form("hfpT_R0%.0lf",fR[r]*10);
		mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} (GeV/c)", nptbins, ptminbin, ptmaxbin)); 
           	hname = Form("hfeta_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet eta;#eta", netabins, etaminbin, etamaxbin));
		hname = Form("hfphi_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "full jet phi;#phi", nphibins, phiminbin, phimaxbin));
		hname = Form("hNF_R0%.0lf",fR[r]*10);
	        mOutList->Add(new TH1D(hname, "jet neutral energy fraction; NEF", 100, 0, 1));
	
	        for (int centbin = 1; centbin < 8; centbin++) {
			//hname = Form("hjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		//mOutList->Add(new TH1D(hname, "jet p_{T}; p_{T} (GeV/c)", nptbins, 0, ptmaxbin));
			//full jet histos
             /*   hname = Form("hDeltaR_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "deltaR vs reco pT",nptbins, ptminbin, ptmaxbin,100, 0, 1));
                hname = Form("hDeltaRw_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "deltaR vs reco pT",nptbins, ptminbin, ptmaxbin,100, 0, 1));
                hname = Form("hNconst_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "Nconst vs reco pT",nptbins, ptminbin, ptmaxbin,150, 0, 150));
                hname = Form("hNconstw_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "Nconst vs reco pT",nptbins, ptminbin, ptmaxbin,150, 0, 150));
                hname = Form("hMCNconst_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "Nconst vs true pT",nptbins, ptminbin, ptmaxbin,150, 0, 150));
                hname = Form("hMCNconstw_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH2D(hname, "Nconst vs true pT",nptbins, ptminbin, ptmaxbin,150, 0, 150));
                hname = Form("hConstpTMCCH_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "MC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTRCCH_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "RC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTMatchedMCCH_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "matched MC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTMatchedRCCH_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "matched RC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTMCN_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "MC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTRCN_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "RC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTMatchedMCN_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "matched MC constituents pT", nptbins, ptminbin, ptmaxbin));
                hname = Form("hConstpTMatchedRCN_R0%.0lf_centbin%i",fR[r]*10, centbin);
                mOutList->Add(new TH1D(hname, "matched RC constituents pT", nptbins, ptminbin, ptmaxbin));*/

            		hname = Form("hfjetpT_R0%.0lf_centbin%i",fR[r]*10, centbin);
		        mOutList->Add(new TH1D(hname, "full jet p_{T}; p_{T} (GeV/c)", nptbins, 0, ptmaxbin));
			hname = Form("hjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
		      //  mOutList->Add(new TH1D(hname, "jet p_{T}^{lead}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
		        hname = Form("hfjetpTlead_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
			hname = Form("hfjetpTleadNeutral_R0%.0lf_centbin%i",fR[r]*10, centbin);
		        mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, N}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
 			hname = Form("hfjetpTleadCharged_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "full jet p_{T}^{lead, Ch}; p_{T} (GeV/c)", 2*npTleadbins, pTleadmin, pTleadmax));
                        hname = Form("hfNjets_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH1D(hname, "number of full jets; N", 100, 0, 100));
            		hname = Form("hetaphi_MCmatched_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, ";eta^{MC}_{jet}; phi^{MC}_{jet}",netabins, etaminbin, etamaxbin,nphibins,phiminbin,phimaxbin));
                       	hname = Form("hetaphi_RCmatched_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, ";eta_{jet}; phi_{jet}",netabins, etaminbin, etamaxbin,nphibins,phiminbin,phimaxbin)); 		
            		hname = Form("hpTleads_R0%.0lf_centbin%i",fR[r]*10, centbin);
            		mOutList->Add(new TH2D(hname, "; p^{det}_{T,lead} (GeV/c); p^{true}_{T,lead} (GeV/c)",2*npTleadbins, pTleadmin, pTleadmax,2*npTleadbins, pTleadmin, pTleadmax)); 		
            		
            		for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                		//hname = Form("hpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
                		TString hdesc = Form("jet p_{T} for p_{T}lead>%i ; p_{T}^{reco} (GeV/c)",pTl);
                		//mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
                		
				hname = Form("hfpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
               			hdesc = Form("full jet p_{T} for p_{T}lead>%i ; p_{T}^{reco} (GeV/c)",pTl);
            			mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
                
              			hname = Form("hMCpT_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
              			hdesc = Form("MC jet p_{T} for p_{T}lead>%i ; p^{MC}_{T} (GeV/c)",pTl);
              			mOutList->Add(new TH1D(hname, hdesc, nptbins, ptminbin, ptmaxbin));
              			
              			hname = Form("hResponseMatrix_pTl%i_R0%.0lf_centbin%i",pTl,fR[r]*10,centbin);
              			hdesc = "; p^{det}_{T} (GeV/c); p^{true} (GeV/c)";
              			mOutList->Add(new TH2D(hname, hdesc, 360, -20, 100, 360, -20, 100));

              			hname = Form("hMCmatchedpT_MCpTl%i_R0%.0lf_centbin%d",pTl,fR[r]*10, centbin);
              			hdesc = "MC matched jets (p_{T,lead} on MC only); p^{true}_{T} (GeV/c)";
              			mOutList->Add(new TH1D(hname, hdesc, 600, 0, ptmaxbin));
              			
              			hname = Form("hMCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[r]*10, centbin);
              			hdesc = "MC matched jets; p^{true}_{T} (GeV/c)";
              			mOutList->Add(new TH1D(hname, hdesc, 600, 0, ptmaxbin));
              			
              			hname = Form("hRCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[r]*10, centbin);
              			hdesc = "RC matched jets; p^{det}_{T} (GeV/c)";
              			mOutList->Add(new TH1D(hname, hdesc, 800, -20, ptmaxbin));
            		}
               	}
							
        } // if (!isMcMode()) {
	} // for(unsigned int r = 0; r < fR.size(); r++) 
		

    if (isMakerMode() == StPicoJetMaker::kWrite) {
        //TODO: Fill trees with required variables
    }

   // mRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
   // mRefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");
  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt="") {
  
  return;
}

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  return kStOK;
}

// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {





	vector<PseudoJet> jetTracks;
	vector<PseudoJet> jetTracks_emb; //tmp for embedding
	vector<PseudoJet> neutraljetTracks; //from bemc towers only
	vector<PseudoJet> fullTracks;
	vector<PseudoJet> MCjetTracks;
	//vector<PseudoJet> MCjetTowers;	

	fRunNumber = mPicoDst->event()->runId();
	int eventId = mPicoDst->event()->eventId(); //eventID
	//cout << "EVENT: "<< eventId << endl;
	int refMult = mPicoDst->event()->refMult();
	//get centrality						
	double vz = mPrimVtx.z();
//	double vr = TMath::Sqrt(mPrimVtx.x()*mPrimVtx.x() + mPrimVtx.y()*mPrimVtx.y());
//	if (vr > 2) cout << "Vr: " << vr << endl; 
	mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(), vz);
	int centrality = mRefmultCorrUtil->centrality9(); //0 = 0-5 %,..., 8 = 70-80 %
    if (centrality==-1) return kStOk; //no centrality
	float Weight = 1.0;
	Weight = mRefmultCorrUtil->weight();
	float weight = Weight*fWeight; //centrality weight * cross section weight
	static_cast<TH1D*>(mOutList->FindObject("hweight"))->Fill(weight);



    if (centrality == 0) centrality = 1; // merge 0-5% and 5-10% into 0-10%
    if (centrality == 8) centrality = 7; // merge 60-70% and 70-80% into 60-80%


    static_cast<TH1D*>(mOutList->FindObject("hcent"))->Fill(centrality, weight);


    //if (centrality > 1) return kStOk; //REMEMBER NOW ONLY CENTRAL
	
	//THIS FUNCTION WILL NOT WORK ON EMBEDDING, UNLESS THE EMBEDDING IS INTO HT EVENTS				
	//if (!FindTriggerTowers(2)) return kStOk; //2 = HT2, don't continue if there is no HT2-trigger tower with sufficient energy
	
	//MC tracks
	int noMCtracks = mPicoDst->numberOfMcTracks();
	for (int i = 0; i < noMCtracks; i++){
		StPicoMcTrack *mctrk = (StPicoMcTrack*)mPicoDst->mcTrack(i);
		if(mctrk->idVtxStart() > 1) continue; //only primary tracks
		int geantId = mctrk->geantId();
		double mcpt = mctrk->pt();
		double mceta = mctrk->eta();
		if ((geantId > 3 && geantId < 7) || fabs(mceta) > 1.0 || mcpt < 0.2) continue;
		double mcp = mctrk->ptot(); 
		TVector3 mcmom = mctrk->p();
   		double mcphi =mcmom.Phi();
  		if(mcphi<0.0) mcphi += 2.0*TMath::Pi(); 
    		if(mcphi>2.0*TMath::Pi()) mcphi -= 2.0*TMath::Pi();

		static_cast<TH2D*>(mOutList->FindObject("hMcEtaPhi"))->Fill(mceta, mcphi, weight); 

		double mcpx,mcpy,mcpz;
		mcpx = mcmom.x();
		mcpy = mcmom.y();
		mcpz = mcmom.z();

		double mcE = mctrk->energy();
		
		PseudoJet inputMcParticle(mcpx, mcpy, mcpz, mcE);
		//PseudoJet inputNeutralMcParticle(mcpx, mcpy, mcpz, mcp); //assume m = 0
		//cout << inputNeutralMcParticle.perp() << endl; 
		if (mctrk->charge() == 0) {inputMcParticle.set_user_index(0);} 
		else inputMcParticle.set_user_index(i);
		MCjetTracks.push_back(inputMcParticle);
	}
	


	//RC part
	
	
	GetCaloTrackMomentum(mPicoDst,mPrimVtx); //fill array Sump with momenta of tracks which are matched to BEMC

    StEmcPosition* mEmcPosition;
    mEmcPosition = new StEmcPosition();

	double TOWE = 0;
	for (int iTow = 0; iTow < 4800; iTow++){ //get btow info
		StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
		if (!towHit || towHit->isBad()) continue; //if the tower is marked as bad or missing info
		int realtowID = towHit->numericIndex2SoftId(iTow);
        if (BadTowerMap[realtowID]) continue; //exclude bad towers (map in JetInfo.h)

      //  cout << "ADC: "<< (towHit->adc()>>4) << endl;

        double towE = GetTowerCalibEnergy(iTow+1); //get tower energy
		TOWE=towE; //just keep track of the original energy for trigger approximation

        if(towErrPlus == true){
            towE = towE + 0.038*towE;
        }
        if(towErrMinus == true){
            towE = towE - 0.038*towE;
        }

		towE-= fHadronCorr*Sump[iTow]; //subtract hadronic energy deposition
		if (towE < 0) towE = 0;
						
		StEmcGeom* mEmcGeom;
		mEmcGeom = StEmcGeom::getEmcGeom("bemc");
		float Toweta_tmp = 0, Towphi = 0;
		mEmcGeom->getEtaPhi(realtowID,Toweta_tmp,Towphi);
        StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(StThreeVectorF(mPrimVtx.x(),mPrimVtx.y(),mPrimVtx.z()), realtowID);
    //    float Toweta2 = vertexCorrectedEta(Toweta_tmp, vz); //max eta 1.05258 max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
        float Toweta = towerPosition.pseudoRapidity();
        static_cast<TH2D*>(mOutList->FindObject("heta_phi_tow"))->Fill(Toweta, Towphi+TMath::Pi(), weight);
		double ET = towE/cosh(Toweta);
		static_cast<TH1D*>(mOutList->FindObject("hET_tow"))->Fill(ET, weight);
		if (ET > 30) {/*cout << towE << endl;*/ continue;} //ignore E > 30 GeV towers 
		//no clustering
		double px,py,pz;
		//px = towE*cos(Towphi)/cosh(Toweta);
		//py = towE*sin(Towphi)/cosh(Toweta);
		px = ET*cos(Towphi);
		py = ET*sin(Towphi);
		pz = towE*tanh(Toweta);
	
		PseudoJet inputTower(px, py, pz, towE);
		if (inputTower.perp() > fETmincut){
		inputTower.set_user_index(0); //default index is -1, 0 means neutral particle
		//THIS LINE WILL NOT WORK
		//if (find(Triggers.begin(), Triggers.end(), realtowID)!=Triggers.end()) inputTower.set_user_index(2); //mark trigger towers with user_index 2
        int ADC = towHit->adc()>>4;
        static_cast<TH1D*>(mOutList->FindObject("hADC"))->Fill(ADC, weight);
		if (ADC > fTrgthresh) inputTower.set_user_index(9999); //mark trigger towers with user_index 9999
		neutraljetTracks.push_back(inputTower);}
	} //end get btow info

    TRandom3 randGen;

	//loop over primary tracks
	for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
        	StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);
        if(trackErr == true) {
            double randomNumber = randGen.Rndm();
            if (randomNumber > 0.96) { continue; }
        }
		double pT = trk->pMom().Perp(); //using primary tracks
		if(pT != pT) continue; // NaN test. 		
        	float eta = trk->pMom().PseudoRapidity();
        	float phi = trk->pMom().Phi();
        	float dca = (mPrimVtx - trk->origin()).Mag();
		float charged = trk->charge();
		
        	static_cast<TH1D*>(mOutList->FindObject("hpT_tr"))->Fill(pT, weight);
        	static_cast<TH2D*>(mOutList->FindObject("heta_phi_tr"))->Fill(phi + TMath::Pi(), eta,  weight);
		static_cast<TH1D*>(mOutList->FindObject("heta_tr"))->Fill(eta, weight);
		static_cast<TH1D*>(mOutList->FindObject("hphi_tr"))->Fill(phi + TMath::Pi(), weight); //to shift by pi
	        static_cast<TH2D*>(mOutList->FindObject("hdca_z_tr"))->Fill(dca, vz, weight);
	        static_cast<TH2D*>(mOutList->FindObject("hdca_pT"))->Fill(dca, pT, weight);
	        static_cast<TH1D*>(mOutList->FindObject("hdca_tr"))->Fill(dca, weight);
		static_cast<TH1D*>(mOutList->FindObject("hcharged_tr"))->Fill(charged, weight);


        //PseudoJet inputParticle(trk->gMom().x(), trk->gMom().y(), trk->gMom().z(), trk->gMom().Mag());
				//primary tracks        
		PseudoJet inputParticle(trk->pMom().x(), trk->pMom().y(), trk->pMom().z(), trk->pMom().Mag());
		if (trk->qaTruth() > 95) {inputParticle.set_user_index(trk->idTruth()-1);/*cout << "MC track pT " << mPicoDst->mcTrack(trk->idTruth()-1)->pt() << " matched to " << pT << endl;*/}
        	jetTracks.push_back(inputParticle);
		} 		//end loop over primary tracks

	//match tracks
	//MatchTracks(MCjetTracks,jetTracks);

	fullTracks = neutraljetTracks;
	fullTracks.insert(fullTracks.end(), jetTracks.begin(), jetTracks.end()); //commenting this line will cause only neutral jets, MAX NEUTRAL FRACTION HAS TO BE TURNED OFF

	//make charged jets
	//background estimation - charged jets
	/*JetDefinition jet_def_bkgd(kt_algorithm, fRBg);
	AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
	if (centrality == 0 || centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
	Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(nJetsRemove)) * SelectorPtMin(0.01);
	JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	bkgd_estimator.set_particles(jetTracks);

	float rho = bkgd_estimator.rho();		
	float rho_sigma = bkgd_estimator.sigma();
	static_cast<TH1D*>(mOutList->FindObject("hrho"))->Fill(rho, weight);
	static_cast<TH2D*>(mOutList->FindObject("hrho_mult"))->Fill(refMult, rho, weight);
	*/
	// full jets
	JetDefinition fjet_def_bkgd(kt_algorithm, fRBg);
	AreaDefinition farea_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
	if (centrality == 1) nJetsRemove = 2;//remove two hardest jets in central collisions, one in others
	Selector fselector =  (!SelectorNHardest(nJetsRemove)) * SelectorAbsEtaMax(1.0) * SelectorPtMin(0.01);
	JetMedianBackgroundEstimator fbkgd_estimator(fselector, fjet_def_bkgd, farea_def_bkgd);
	fbkgd_estimator.set_particles(fullTracks);

	float frho   = fbkgd_estimator.rho();
	//float frho_sigma = fbkgd_estimator.sigma();
	static_cast<TH1D*>(mOutList->FindObject("hfrho"))->Fill(frho, weight);
	static_cast<TH2D*>(mOutList->FindObject("hfrho_mult"))->Fill(refMult, frho, weight);


    	for (unsigned int i = 0; i < fR.size(); i++) { 
        	float maxRapJet = mPicoCuts->getCutEta() - fR[i];

        	if (isMcMode()) {	//not used

        	} else { //!isMcMode()
	        float maxRapJet = 1 - fR[i];
		//MC jets first
		
		//setup fastjet
		JetDefinition Mcjet_def(antikt_algorithm, fR[i]);
        	AreaDefinition Mcarea_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
            //run jet reconstruction
		ClusterSequenceArea Mcclust_seq_hard(MCjetTracks, Mcjet_def, Mcarea_def);
		vector<PseudoJet> Mcjets_all = sorted_by_pt(Mcclust_seq_hard.inclusive_jets(fJetPtMin));
          	Selector McFiducial_cut_selector = SelectorAbsEtaMax(maxRapJet)* SelectorPtMin(0.01)* SelectorPtMax(1.5*fpThatmax); //throw out jets with pT larger than 3*pThat to eliminate high-weight fluctuations; // Fiducial cut for jets
            	vector<PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
		vector<double> McPtLeads;
		//int nMCjets10 = 0;


        int NumberOfConst = 0;
		for(unsigned int pjet = 0; pjet < McJets.size(); pjet++) {
                	float phi_jet = McJets[pjet].phi();
                	float eta_jet = McJets[pjet].eta();
               		float pT_jet = McJets[pjet].perp();
                	float area_jet = McJets[pjet].area();
               	 	vector<PseudoJet> constituents = sorted_by_pt(McJets[pjet].constituents());							
                	float pTlead = constituents[0].perp();
            for(unsigned int icc = 0; icc < constituents.size(); ++icc) {
                if (constituents[icc].perp()>0.2) {
                   /* if(pT_jet>20.0){
                        if (constituents[icc].user_index() == 0 || constituents[icc].user_index() == 9999) {

                            static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMCN_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituents[icc].perp());
                            }
                        else {
                            static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMCCH_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituents[icc].perp());
                        }
                        }*/
                    NumberOfConst++;
                }
            }


                McPtLeads.push_back(pTlead);
                	for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                        	if(pTl < pTlead) static_cast<TH1D*>(mOutList->FindObject(Form("hMCpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_jet, weight);
                    	}
                    	//if (pT_jet > 10) nMCjets10++;
  
            	} // for(unsigned int pjet = 0; pjet < McJets.size(); pjet++)
        	
        	
        	//RC jets
		//setup fastjet
		//full jet reconstruction
		JetDefinition fjet_def(antikt_algorithm, fR[i]);
		AreaDefinition farea_def(active_area_explicit_ghosts, GhostedAreaSpec(fGhostMaxrap, 1, 0.01));
        	//run full jet reconstruction
		ClusterSequenceArea fclust_seq_hard(fullTracks, fjet_def, farea_def);
            	vector<PseudoJet> fjets_all = sorted_by_pt(fclust_seq_hard.inclusive_jets(fJetPtMin));
                Selector fFiducial_cut_selector = SelectorAbsEtaMax(maxRapJet); //throw out jets with pT larger than 1.5*pThat (upper edge) to eliminate high-weight fluctuations;// Fiducial cut for jets
            	vector<PseudoJet> fjets = fFiducial_cut_selector(fjets_all);

		int naccJets = 0;
		vector<double> RcPtLeads;
		vector<PseudoJet> RcJets;
		for(unsigned int pjet = 0; pjet < fjets.size(); pjet++) {
		        vector<PseudoJet> constituents = sorted_by_pt(fjets[pjet].constituents());
			//look only for jets with associated trigger tower
			bool istriggerjet = false;
			for(unsigned int ic = 0; ic < constituents.size(); ++ic) {
				if (constituents[ic].user_index() == 9999) {istriggerjet = true; break;}
			}
			if (!istriggerjet) continue;
                	float phi_jet = fjets[pjet].phi();
                	float eta_jet = fjets[pjet].eta();
                	float pT_jet = fjets[pjet].perp();
                	float area_jet = fjets[pjet].area();
                	int nparticles = constituents.size();
			//totaljetE += fjets[pjet].E();
			float pTcorr_jet = pT_jet - area_jet*frho;
			float pTlead = constituents[0].perp();
			
			RcPtLeads.push_back(pTlead);
	                RcJets.push_back(fjets[pjet]);
						
                	//set acceptance
                	float etaMinCut = -(maxRapJet);
                	float etaMaxCut = (maxRapJet);

                	if(eta_jet < etaMinCut || eta_jet > etaMaxCut) continue; // fiducial acceptance

                	static_cast<TH1D*>(mOutList->FindObject(Form("hfjetarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, weight);
              		static_cast<TH2D*>(mOutList->FindObject(Form("hfjetpTarea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pT_jet, weight);
                	static_cast<TH2D*>(mOutList->FindObject(Form("hfjetpTcorrArea_R0%.0lf",fR[i]*10)))->Fill(area_jet, pTcorr_jet, weight);

                	if(area_jet < fAcuts[i]) continue;
			
			float neutralpT = 0;
			//is the leading particle charged or neutral?
			for(unsigned int ic = 0; ic < constituents.size(); ++ic) {
				if (constituents[ic].user_index() == 0 || constituents[ic].user_index()== 9999) neutralpT+=constituents[ic].pt(); //select towers
				float ceta = constituents[ic].eta();
				float cphi = constituents[ic].phi();
				static_cast<TH1D*>(mOutList->FindObject(Form("hfceta_R0%.0lf",fR[i]*10)))->Fill(ceta, weight);
				static_cast<TH1D*>(mOutList->FindObject(Form("hfcphi_R0%.0lf",fR[i]*10)))->Fill(cphi, weight);
                if (constituents[ic].perp()>0.2) {
                   /* if(pT_jet>20.0){
                        if (constituents[ic].user_index() == 0 || constituents[ic].user_index() == 9999) {

                            static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTRCN_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituents[ic].perp());
                        }
                        else {
                            static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTRCCH_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituents[ic].perp());
                        }
                    }*/
                    NumberOfConst++;
                }
			}
        //    cout << "Number of constituents "<< NumberOfConst <<endl;




			float nfraction = neutralpT/pT_jet;
			//cout << "neutral fraction " << nfraction << endl;
			static_cast<TH1D*>(mOutList->FindObject(Form("hNF_R0%.0lf",fR[i]*10)))->Fill(nfraction, weight);
			if (nfraction > maxneutralfrac) continue; //keep only jets with reasonable neutral energy fraction of jet pT (default 0.95)	
			naccJets++;		
							
			float leadindex = constituents[0].user_index();
			//cout << leadindex << endl;
			if (leadindex == 0 || leadindex == 9999) { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadNeutral_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);} 
			else { static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTleadCharged_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);}
			static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpT_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pT_jet, weight);
                	static_cast<TH1D*>(mOutList->FindObject(Form("hfjetpTlead_R0%.0lf_centbin%d",fR[i]*10, centrality)))->Fill(pTlead, weight);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_R0%.0lf",fR[i]*10)))->Fill(pT_jet, weight);
			//static_cast<TH2D*>(mOutList->FindObject(Form("heta_phi_R0%.0lf", fR[i]*10)))->Fill(eta_jet, phi_jet);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfeta_R0%.0lf", fR[i]*10)))->Fill(eta_jet, weight);
			static_cast<TH1D*>(mOutList->FindObject(Form("hfphi_R0%.0lf", fR[i]*10)))->Fill(phi_jet, weight);
	                static_cast<TH2D*>(mOutList->FindObject(Form("hfnparticlesinjet_R0%.0lf",fR[i]*10)))->Fill(nparticles, pTlead);

	                for(Int_t pTl = 0; pTl < npTlead; pTl++) {
                     	   if(pTl < pTlead) {static_cast<TH1D*>(mOutList->FindObject(Form("hfpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pTcorr_jet, weight);}
                    	}


		} // for(unsigned int pjet = 0; pjet < jets.size(); pjet++)
            	static_cast<TH1D*>(mOutList->FindObject(Form("hfNjets_R0%.0lf_centbin%d", fR[i]*10,centrality)))->Fill(naccJets, weight);
           	//cout << "total jet energy in this event: " << totaljetE << endl;
   		if (RcJets.size() == 0) continue;
		vector<pair<PseudoJet, PseudoJet>> Matched;
		vector<pair<double, double>> MatchedpTleads;
		vector<pair<double, double>> MatchedNeutralFraction;
        vector<double> differEta, differPhi;
        vector<pair<double, double>> deltaR;
		//vector<pair<int, int>> MatchedNNeutral, MatchedNCharged, MatchedNTot;
		MatchJetsEtaPhi(McJets, RcJets, McPtLeads, RcPtLeads, &Matched, &MatchedpTleads, &MatchedNeutralFraction, /*&MatchedNNeutral, &MatchedNCharged, &MatchedNTot, */fR[i], &differEta, &differPhi, &deltaR);
		//cout << deltaR << " " << deltapT << " " << pTtrue << endl;
                for (double value : differPhi) {
                    static_cast<TH1D*>(mOutList->FindObject(Form("hphi_MCRC_R0%.0lf", fR[i]*10)))->Fill(value + TMath::Pi());
                    static_cast<TH1D*>(mOutList->FindObject(Form("hphi_MCRCw_R0%.0lf", fR[i]*10)))->Fill(value + TMath::Pi(), weight);
                }

                for (double value : differEta) {
                    static_cast<TH1D*>(mOutList->FindObject(Form("heta_MCRC_R0%.0lf", fR[i]*10)))->Fill(value);
                    static_cast<TH1D*>(mOutList->FindObject(Form("heta_MCRCw_R0%.0lf", fR[i]*10)))->Fill(value, weight);

                }

                for (size_t j = 0; j < differEta.size(); ++j) {
                    double etaValue = differEta[j];
                    double phiValue = differPhi[j];

                    static_cast<TH2D*>(mOutList->FindObject(Form("hEtaPhi_MC-RC_R0%.0lf", fR[i]*10)))->Fill(etaValue, phiValue + TMath::Pi());
                    static_cast<TH2D*>(mOutList->FindObject(Form("hEtaPhi_MC-RCw_R0%.0lf", fR[i]*10)))->Fill(etaValue, phiValue + TMath::Pi(), weight);
                }

                for (int j = 0; j < deltaR.size(); ++j) {
                    double deltaRvalue = deltaR[j].first;
                    double pTvalue = deltaR[j].second;

//                    static_cast<TH2D*>(mOutList->FindObject(Form("hDeltaR_R0%.0lf", fR[i]*10)))->Fill(pTvalue, deltaRvalue);
                  //  static_cast<TH2D*>(mOutList->FindObject(Form("hDeltaR_R0%.0lf_centbin%i", fR[i]*10, centrality)))->Fill(pTvalue, deltaRvalue);
                  //  static_cast<TH2D*>(mOutList->FindObject(Form("hDeltaRw_R0%.0lf_centbin%i", fR[i]*10, centrality)))->Fill(pTvalue, deltaRvalue, weight);

                }
                for (unsigned int j = 0; j < Matched.size(); j++) {
			double pT_det = Matched[j].second.perp();
			double pT_true = Matched[j].first.perp();
			double pT_corr_det = pT_det - Matched[j].second.area()*frho;
			double MCmatchedeta = Matched[j].first.eta();
			double MCmatchedphi = Matched[j].first.phi();
			double RCmatchedeta = Matched[j].second.eta();
			double RCmatchedphi = Matched[j].second.phi();

            int NumberOfConstMC = 0;
            int NumberOfConstRC = 0;

                vector<PseudoJet> constituentsMC = sorted_by_pt(Matched[j].first.constituents());
                    for(unsigned int ic = 0; ic < constituentsMC.size(); ++ic) {
                        if (constituentsMC[ic].perp()>0.2) {
                           /* if(pT_true>20.0) {
                                if (constituentsMC[ic].user_index() == 0 || constituentsMC[ic].user_index() == 9999) {

                                    static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMatchedMCN_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituentsMC[ic].perp());
                                }
                                else {
                                    static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMatchedMCCH_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituentsMC[ic].perp());
                                }
                            }*/
                            NumberOfConstMC++;
                        }
                    }
                vector<PseudoJet> constituentsRC = sorted_by_pt(Matched[j].second.constituents());
                    for(unsigned int icc = 0; icc < constituentsRC.size(); ++icc) {
                        if (constituentsRC[icc].perp()>0.2) {
                           /* if(pT_corr_det>20.0) {
                                if (constituentsRC[icc].user_index() == 0 || constituentsRC[icc].user_index() == 9999) {

                                    static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMatchedRCN_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituentsRC[icc].perp());
                                }
                                else {
                                    static_cast<TH1D *>(mOutList->FindObject(Form("hConstpTMatchedRCCH_R0%.0lf_centbin%i", fR[i] * 10, centrality)))->Fill(constituentsRC[icc].perp());
                                }
                            }*/
                            NumberOfConstRC++;
                        }
                    }


        /*    static_cast<TH2D*>(mOutList->FindObject(Form("hMCNconst_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(pT_true, NumberOfConstMC);
            static_cast<TH2D*>(mOutList->FindObject(Form("hMCNconstw_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(pT_true, NumberOfConstMC, weight);

            static_cast<TH2D*>(mOutList->FindObject(Form("hNconst_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(pT_corr_det, NumberOfConstRC);
            static_cast<TH2D*>(mOutList->FindObject(Form("hNconstw_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(pT_corr_det, NumberOfConstRC, weight);*/

			static_cast<TH2D*>(mOutList->FindObject(Form("hetaphi_MCmatched_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(MCmatchedeta,MCmatchedphi, weight);
			static_cast<TH2D*>(mOutList->FindObject(Form("hetaphi_RCmatched_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(RCmatchedeta,RCmatchedphi, weight);	

			//double pTlead = MatchedpTleads[j].second; //pTlead cut on detector level only
			double MatchedMCpTlead = MatchedpTleads[j].first;
			double MatchedRCpTlead = MatchedpTleads[j].second;	
			//cout << "MC: " << MatchedMCpTlead << " det: " << MatchedRCpTlead; {MatchedRCpTlead < MatchedMCpTlead ? cout << " RC " << endl : cout << " MC " << endl; }
			//double pTlead = min(MatchedMCpTlead,MatchedRCpTlead); //pTlead cut on both levels 
			double pTlead = MatchedRCpTlead; //pTlead cut on detector level only
			double matchedNF = MatchedNeutralFraction[j].first;

                    if (matchedNF < 0.01) continue; //throw out track-only MC jets
              /*      cout << "Eta MC  " << MCmatchedeta << "Eta RC " << RCmatchedeta << endl;
                    cout << "Phi MC  " << MCmatchedphi << "Phi RC " << RCmatchedphi << endl;
                    cout << "pT  MC  " << pT_true  << "pT RC  " << pT_det <<"pT_corr  "<<pT_corr_det<<endl;
                    cout << "Matched neutral fraction" << matchedNF << endl;*/

                    static_cast<TH2D*>(mOutList->FindObject(Form("hpTleads_R0%.0lf_centbin%i",fR[i]*10, centrality)))->Fill(MatchedRCpTlead,MatchedMCpTlead,weight);
			for(Int_t pTl = 0; pTl < npTlead; pTl++) {
				if(pTl < pTlead) {
				static_cast<TH2D*>(mOutList->FindObject(Form("hResponseMatrix_pTl%i_R0%.0lf_centbin%i",pTl,fR[i]*10,centrality)))->Fill(pT_corr_det, pT_true, weight);
				static_cast<TH1D*>(mOutList->FindObject(Form("hRCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_corr_det, weight);
				static_cast<TH1D*>(mOutList->FindObject(Form("hMCmatchedpT_pTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_true, weight);
				} 
			 	if (pTl < MatchedpTleads[j].first){
                     //ZDE
			 	static_cast<TH1D*>(mOutList->FindObject(Form("hMCmatchedpT_MCpTl%i_R0%.0lf_centbin%d",pTl,fR[i]*10, centrality)))->Fill(pT_true, weight);
			 	}
			}
		}
		Matched.clear();
		MatchedpTleads.clear();
		MatchedNeutralFraction.clear();
       		} //if (!isMcMode())
	} //for (unsigned int i = 0; i < fR.size(); i++)

	//clear vectors/arrays
	Sump.fill(0);
	Triggers.clear();
	
	return kStOK;
}

Double_t StPicoHFJetMaker::GetTowerCalibEnergy(Int_t TowerId)
{
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(TowerId-1));

  Float_t pedestal, rms;
  Int_t status;
  
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);	

  Double_t *TowerCoeff;
  if(fRunNumber <= 15094020) TowerCoeff = CPre;
  else TowerCoeff = CLowMidHigh;
  
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal);

  return calibEnergy;
}



//----------------------------------------------------------------------------- 
////Correct tower eta for Vz position //// Not used anymore
//----------------------------------------------------------------------------- 
Double_t StPicoHFJetMaker::vertexCorrectedEta(double eta, double vz) {
    double tower_theta = 2.0 * atan(exp(-eta));
    double z = 0.0;
    if (eta != 0.0) z = mBarrelRadius / tan(tower_theta);
    double z_diff = z - vz;
    double theta_corr = atan2(mBarrelRadius, z_diff);
    double eta_corr = -log(tan(theta_corr / 2.0));
    return eta_corr;
}

//----------------------------------------------------------------------------- 
//Fill array with momentum of BEMC-matched tracks
//----------------------------------------------------------------------------- 
Bool_t StPicoHFJetMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst, TVector3 mPrimVtx) {
	//loop over global tracks  - towers

	UInt_t nTracks = mPicoDst->numberOfTracks();
	//cout << nTracks << endl;
	for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
        StPicoTrack *trk = mPicoDst->track(itrack);
	TVector3 gMom = trk->gMom();
	//using global tracks
	double pT = gMom.Perp();
	if(pT != pT || pT < 0.2) continue;
        float eta = gMom.PseudoRapidity();
	if (fabs(eta) > 1) continue;
        float phi = gMom.Phi();

	float nHitsFit = trk->nHitsFit();
	float nHitsMax = trk->nHitsMax();
	if (nHitsFit < 15 || nHitsFit/nHitsMax < 0.52) continue; //some basic QA cuts
	double Bfield = mPicoDst->event()->bField();

	StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
	float vtx_x = mPrimVtx.x();
	float vtx_y = mPrimVtx.y();
	float vtx_z = mPrimVtx.z();

        float dca_z = abs(trk->gDCAz(mPicoDst->event()->primaryVertex().z()));
        if (fabs(dca_z) > maxdcazhadroncorr) continue; 
	int TowIndex = -99999;
	TowIndex = trk->bemcTowerIndex();
	float p = 0;
	if (TowIndex >= 0) {
//        if (TowIndex % 20 == 0){
//            cout << "TowIndex " << TowIndex << " Eta " << eta << endl;


//        }

        p = gMom.Mag();
        Sump[TowIndex] += p;
		//cout << p << " " << Sump[TowIndex-1] << endl;
		}
	}// END global track loop
	
	return true;
}

//----------------------------------------------------------------------------- 
//Fill array with ID of towers which are marked as HTlevel and are above the desired trigger threshold 
// IN DIJET EMBEDDING INTO MB, THERE IS NO HT TRIGGERS
//----------------------------------------------------------------------------- 
 Int_t StPicoHFJetMaker::FindTriggerTowers(Int_t level = 2) {
 
 	if (level < 1 || level > 3) {cout << "Wrong trigger level, this function cannot be used" << endl; return -1;}
 	
 	UInt_t ntrg = mPicoDst->numberOfEmcTriggers();
	for (int i = 0; i < ntrg; i++){ //get btow info
		StPicoEmcTrigger *trg = mPicoDst->emcTrigger(i);
		if (!trg) {/*cout << "no trigger info" << endl;*/ continue;} 
		if ((level == 2 && !trg->isHT2()) || (level == 1 && !trg->isHT1()) || (level == 3 && !trg->isHT3())) {/*cout << "not HT2 trigger" << endl;*/ continue;}
		int towid = trg->id();
		if (BadTowerMap[towid]) continue; 
		float energy = GetTowerCalibEnergy(towid);
        float ADC = mPicoDst->btowHit(towid-1)->adc();
		if (ADC > fTrgthresh) Triggers.push_back(towid);
	} 	
	
	return Triggers.size();
 }


