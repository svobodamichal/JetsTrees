#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>

// ROOT 5 / CINT compatible: C++03 syntax, no C++11 smart pointers.
void rawComparisonFatal(const char* message) {
  std::cerr << "Comparison failed: " << message << std::endl;
  gSystem->Exit(1);
}

// Raw = detector-level, background-corrected pT (reco_pt_corr).
// Counts retain the centralityWeight applied by make_hists_data.C.
void compare_raw_calibrations(const char* origFile="hists_data_OrigCalib.root",
                             const char* hsFile="hists_data_HSCalib.root",
                             const char* outDir="raw_calibration_comparison",
                             int rebin=1, double xmin=-20., double xmax=60.) {
  {
    if (rebin<1 || xmin>=xmax) rawComparisonFatal("Invalid rebin or x range");
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);
    TFile* fo = TFile::Open(origFile,"READ");
    TFile* fh = TFile::Open(hsFile,"READ");
    if (!fo || fo->IsZombie() || !fh || fh->IsZombie())
      rawComparisonFatal("Cannot open one of the input ROOT files");
    if (gSystem->mkdir(outDir,kTRUE)!=0 && gSystem->AccessPathName(outDir))
      rawComparisonFatal("Cannot create output directory");
    TString output= TString(outDir)+"/raw_calibration_comparison.root";
    if (output==origFile || output==hsFile) rawComparisonFatal("Output matches input");
    TFile* fout = TFile::Open(output.Data(),"RECREATE");
    if (!fout || fout->IsZombie()) rawComparisonFatal("Cannot create output ROOT file");
    const int leads[]={0,5,7,9};
    int made=0;
    TIter ir(fo->GetListOfKeys()); TKey* kr;
    while ((kr=(TKey*)ir())) {
      TString rn=kr->GetName(); double radius=0;
      TDirectory* dr=fo->GetDirectory(rn);
      if (!dr || sscanf(rn.Data(),"R%lf",&radius)!=1) continue;
      TDirectory* outR=fout->mkdir(rn);
      TIter ic(dr->GetListOfKeys()); TKey* kc;
      while ((kc=(TKey*)ic())) {
        TString cn=kc->GetName(); TDirectory* dc=dr->GetDirectory(cn);
        if (!dc) continue;
        TDirectory* outC=outR->mkdir(cn);
        TString plotDir=TString(outDir)+"/"+rn+"/"+cn;
        gSystem->mkdir(plotDir,kTRUE);
        for (int j=0;j<4;++j) {
          TString hn=Form("hSpec_recoPtCorr_ptlead%d_%s_%s",leads[j],rn.Data(),cn.Data());
          TH1* ao=dynamic_cast<TH1*>(dc->Get(hn));
          TH1* ah=dynamic_cast<TH1*>(fh->Get(rn+"/"+cn+"/"+hn));
          if (!ao || !ah) rawComparisonFatal(Form("Missing histogram: %s/%s/%s",rn.Data(),cn.Data(),hn.Data()));
          if (ao->GetDimension()!=1 || ah->GetDimension()!=1 || ao->GetNbinsX()!=ah->GetNbinsX())
            rawComparisonFatal("Histogram dimensions/binning differ");
          for (int b=1;b<=ao->GetNbinsX()+1;++b)
            if (std::fabs(ao->GetXaxis()->GetBinLowEdge(b)-ah->GetXaxis()->GetBinLowEdge(b))>1e-9)
              rawComparisonFatal("Histogram bin edges differ");
          if (ao->GetNbinsX()%rebin) rawComparisonFatal("Rebin must divide the number of bins (160)");
          TH1* o = (TH1*)ao->Clone(Form("original_ptlead%d",leads[j]));
          TH1* h = (TH1*)ah->Clone(Form("HS_ptlead%d",leads[j]));
          if (rebin>1) {o->Rebin(rebin);h->Rebin(rebin);}
          o->Scale(1.,"width"); h->Scale(1.,"width");
          TGraphErrors* ratio = new TGraphErrors(); ratio->SetName(Form("ratio_Orig_over_HS_ptlead%d",leads[j]));
          double lo=1e300, hi=0, rlo=1, rhi=1;
          for (int b=1;b<=o->GetNbinsX();++b) {
            double x=o->GetBinCenter(b), a=o->GetBinContent(b), d=h->GetBinContent(b);
            if (x<xmin || x>xmax) continue;
            if (a>0) {lo=std::min(lo,a);hi=std::max(hi,a);}
            if (d>0) {lo=std::min(lo,d);hi=std::max(hi,d);}
            if (d<=0) continue; // Undefined ratios are omitted, never shown as zero.
            double ea=o->GetBinError(b), ed=h->GetBinError(b), y=a/d;
            double e=std::sqrt((ea/d)*(ea/d)+(a*ed/(d*d))*(a*ed/(d*d)));
            int n=ratio->GetN();ratio->SetPoint(n,x,y);ratio->SetPointError(n,0,e);
            rlo=std::min(rlo,y-e);rhi=std::max(rhi,y+e);
          }
          if (hi<=0) rawComparisonFatal("No positive spectrum bins in requested range");
          TCanvas* c = new TCanvas(Form("comparison_ptlead%d",leads[j]),"Raw calibration comparison",850,800);
          TPad* top = new TPad("top","",0,.28,1,1);
          TPad* bottom = new TPad("bottom","",0,0,1,.28);
          top->SetLeftMargin(.17);top->SetRightMargin(.04);top->SetBottomMargin(.025);top->SetTopMargin(.05);
          bottom->SetLeftMargin(.17);bottom->SetRightMargin(.04);bottom->SetTopMargin(.025);bottom->SetBottomMargin(.32);
          top->Draw();bottom->Draw();top->cd();top->SetLogy();top->SetTicks(1,1);
          o->SetTitle("");o->SetLineColor(kRed+1);o->SetMarkerColor(kRed+1);o->SetMarkerStyle(24);o->SetMarkerSize(.65);
          h->SetLineColor(kBlue+1);h->SetMarkerColor(kBlue+1);h->SetMarkerStyle(20);h->SetMarkerSize(.65);
          o->GetXaxis()->SetRangeUser(xmin,xmax);o->GetXaxis()->SetLabelSize(0);o->GetXaxis()->SetTitle("");
          o->GetYaxis()->SetTitle("Weighted jet counts / (GeV/#it{c})");
          o->GetYaxis()->SetTitleSize(.045);o->GetYaxis()->SetTitleOffset(1.6);o->GetYaxis()->SetLabelSize(.043);
          o->SetMinimum(lo*.3);o->SetMaximum(hi*20);o->Draw("E1");h->Draw("E1 SAME");
          TLegend* leg = new TLegend(.56,.77,.94,.92);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->SetTextSize(.034);
          leg->AddEntry(o,"Original calibration","lep");leg->AddEntry(h,"Hanseul calibration (baseline)","lep");leg->Draw();
          TString cent=cn;cent.ReplaceAll("CENT_","");cent.ReplaceAll("MID_","");cent.ReplaceAll("PERI_","");cent.ReplaceAll("_","-");
          TLatex tx;tx.SetNDC();tx.SetTextFont(42);tx.SetTextSize(.035);
          tx.DrawLatex(.21,.91,Form("Au+Au %s%%, #sqrt{s_{NN}} = 200 GeV",cent.Data()));
          tx.DrawLatex(.21,.85,Form("anti-#it{k}_{T}, R = %.1f",radius));
          tx.DrawLatex(.21,.79,Form("#it{p}_{T,min}^{lead} = %d GeV/#it{c}",leads[j]));
          bottom->cd();bottom->SetTicks(1,1);
          TH1* frame = (TH1*)o->Clone("ratio_frame");frame->Reset();
          double span=std::max(.4,rhi-rlo);
          frame->SetMinimum(std::min(.5,rlo-.12*span));frame->SetMaximum(std::max(1.5,rhi+.12*span));
          frame->GetYaxis()->SetTitle("OrigCalib / HSCalib");frame->GetYaxis()->SetTitleSize(.105);frame->GetYaxis()->SetTitleOffset(.68);
          frame->GetYaxis()->SetLabelSize(.10);frame->GetYaxis()->SetNdivisions(505);
          frame->GetXaxis()->SetTitle("#it{p}_{T,jet}^{reco,corr} [GeV/#it{c}]");frame->GetXaxis()->SetTitleSize(.12);
          frame->GetXaxis()->SetLabelSize(.105);frame->GetXaxis()->SetTitleOffset(1.05);frame->Draw("AXIS");
          TLine* unity = new TLine(xmin,1,xmax,1);unity->SetLineStyle(2);unity->Draw();
          ratio->SetLineColor(kBlue+1);ratio->SetMarkerColor(kBlue+1);ratio->SetMarkerStyle(20);ratio->SetMarkerSize(.65);ratio->Draw("P SAME");
          TString base=plotDir+Form("/raw_comparison_%s_%s_ptlead%d",rn.Data(),cn.Data(),leads[j]);
          c->SaveAs(base+".pdf");c->SaveAs(base+".png");
          outC->cd();o->Write();h->Write();ratio->Write();c->Write();++made;
          // Canvas owns heap-allocated pads. Pads own all drawn objects below.
          // Write everything first, then destroy the complete plot once.
          o->SetBit(TObject::kCanDelete); h->SetBit(TObject::kCanDelete);
          frame->SetBit(TObject::kCanDelete); ratio->SetBit(TObject::kCanDelete);
          leg->SetBit(TObject::kCanDelete); unity->SetBit(TObject::kCanDelete);
          delete c;
        }
      }
    }
    if (!made) rawComparisonFatal("No matching radius/centrality directories found");
    fout->Close(); fo->Close(); fh->Close();
    delete fout; delete fo; delete fh;
    std::cout<<"Created "<<made<<" comparisons in "<<outDir<<std::endl;
  }
}
