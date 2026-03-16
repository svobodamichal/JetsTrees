// make_efficiency_ratio_pdfs.C  (ROOT5 / old C++ / CINT-friendly)
//
// Creates ratio PDFs: (pThat bin)/(reference pThat bin) for three efficiencies.
// - match: h_match_eff_truth  vs pT^MC
// - trig : h_trig_eff_reco    vs pT^reco,corr
// - pur  : h_pur_eff_reco     vs pT^reco,corr
//
// Assumes your efficiency ROOT files contain top-level tag dirs:
//   R0.2_CENT_0_10_ptlead0, R0.2_CENT_0_10_ptlead5, ...
//
// Output: outdir/{match,trig,pur}/{TAG}/ratio_<type>__<file>__over__<ref>.pdf
//
// Usage:
//   root -l -b -q 'make_efficiency_ratio_pdfs.C("analysis/efficiencies","PDFs_Ratios","50_-1")'

#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TList.h"
#include "TCollection.h"  
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRegexp.h"

#include <vector>
#include <string>
#include <iostream>

using std::cout;
using std::endl;

static void ApplyStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);

  gStyle->SetTitleSize(0.045, "XY");
  gStyle->SetLabelSize(0.040, "XY");
  gStyle->SetTitleOffset(1.10, "X");
  gStyle->SetTitleOffset(1.20, "Y");
}

static TString BaseNoExt(const TString& path)
{
  TString b = gSystem->BaseName(path);
  if (b.EndsWith(".root")) b.ReplaceAll(".root", "");
  return b;
}

static void MkdirP(const TString& d)
{
  if (d.Length() == 0) return;
  gSystem->mkdir(d, kTRUE);
}

// Fetch histogram (returns a heap clone with SetDirectory(0), or 0)
static TH1D* GetH(TFile* f, const TString& tag, const TString& hname)
{
  if (!f) return 0;
  TDirectory* d = (TDirectory*)f->Get(tag);
  if (!d) return 0;

  TH1* h = (TH1*)d->Get(hname);
  if (!h) return 0;

  TH1D* hc = (TH1D*)h->Clone(Form("%s__tmp", hname.Data()));
  hc->SetDirectory(0);
  return hc;
}

// Build ratio safely bin-by-bin with Gaussian propagation
static TH1D* MakeRatio(const TH1D* a, const TH1D* b, const TString& name)
{
  if (!a || !b) return 0;
  if (a->GetNbinsX() != b->GetNbinsX()) return 0;

  TH1D* r = (TH1D*)a->Clone(name);
  r->SetDirectory(0);
  r->Reset("ICES");

  for (int i = 1; i <= a->GetNbinsX(); ++i) {
    double A  = a->GetBinContent(i);
    double eA = a->GetBinError(i);
    double B  = b->GetBinContent(i);
    double eB = b->GetBinError(i);

    if (B <= 0 || A < 0) {
      r->SetBinContent(i, 0.0);
      r->SetBinError(i, 0.0);
      continue;
    }

    double R = A / B;

    double rel2 = 0.0;
    if (A > 0) rel2 += (eA*eA)/(A*A);
    if (B > 0) rel2 += (eB*eB)/(B*B);
    double eR = R * std::sqrt(rel2);

    r->SetBinContent(i, R);
    r->SetBinError(i, eR);
  }

  return r;
}

static void DrawSave(TH1D* h,
                     const TString& outPdf,
                     const TString& xTitle,
                     const TString& header1,
                     const TString& header2)
{
  if (!h) return;

  MkdirP(gSystem->DirName(outPdf));

  TCanvas* c = new TCanvas("c", "c", 900, 700);
  c->SetTicks(1,1);

  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.0);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle("ratio to reference");
  h->GetXaxis()->SetRangeUser(-40, 60);

  // y-range: auto but keep readable around 1
  double ymin =  1e9, ymax = -1e9;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double v = h->GetBinContent(i);
    double e = h->GetBinError(i);
    if (v <= 0) continue;
    if (v - e < ymin) ymin = v - e;
    if (v + e > ymax) ymax = v + e;
  }
  if (ymin < 1e8 && ymax > -1e8) {
    double pad = 0.15*(ymax - ymin);
    if (pad <= 0) pad = 0.10;
    ymin -= pad; ymax += pad;
    if (fabs(ymax - ymin) < 0.05) { ymin = 0.95; ymax = 1.05; }
    if (ymin < 0.0) ymin = 0.0;
    if (ymax > 2.0) ymax = 2.0;
    h->GetYaxis()->SetRangeUser(ymin, ymax);
  } else {
    h->GetYaxis()->SetRangeUser(0.8, 1.2);
  }

  h->Draw("E1");

  TLine* l1 = new TLine(h->GetXaxis()->GetXmin(), 1.0,
                        h->GetXaxis()->GetXmax(), 1.0);
  l1->SetLineStyle(2);
  l1->SetLineWidth(2);
  l1->Draw("same");

  TLatex lat;
  lat.SetNDC(kTRUE);
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.14, 0.94, header1);
  lat.DrawLatex(0.14, 0.90, header2);

  c->SaveAs(outPdf);

  delete l1;
  delete c;
}

void make_efficiency_ratio_pdfs(const char* indir = "analysis/efficiencies",
                                const char* outdir = "analysis/efficiencies/PDFs_Ratios",
                                const char* refPattern = "50_-1")
{
  ApplyStyle();

  TString inDir(indir);
  TString outDir(outdir);
  TString refPat(refPattern);

// collect files (ROOT5/CINT-safe: gSystem->OpenDirectory)
std::vector<TString> files;
TString refFile = "";

void* dirp = gSystem->OpenDirectory(inDir);
if (!dirp) {
  cout << "ERROR: cannot open directory " << inDir << endl;
  return;
}

const char* entry = 0;
while ((entry = gSystem->GetDirEntry(dirp))) {
  TString name(entry);

  // skip . and ..
  if (name == "." || name == "..") continue;

  // only root files
  if (!name.EndsWith(".root")) continue;

  // only your efficiency files
  if (!name.BeginsWith("efficiencies_")) continue;

  // skip full (optional)
  if (name.Contains("full")) continue;

  TString full = inDir + "/" + name;
  files.push_back(full);

  if (name.Contains(refPat)) refFile = full;
}
gSystem->FreeDirectory(dirp);

if (files.empty()) {
  cout << "ERROR: no efficiencies_*.root found in " << inDir << endl;
  return;
}
if (refFile.Length() == 0) {
  cout << "ERROR: reference not found (pattern '" << refPat << "')" << endl;
  for (size_t i=0;i<files.size();++i) cout << "  " << files[i] << endl;
  return;
}

  cout << "Reference: " << refFile << endl;

  // Open reference once
  TFile* fRef = TFile::Open(refFile, "READ");
  if (!fRef || fRef->IsZombie()) { cout << "ERROR: cannot open ref\n"; return; }

  // Define efficiency types
  const int NTYPE = 3;
  TString typeKey[NTYPE]  = {"match","trig","pur"};
  TString hname[NTYPE]    = {"h_match_eff_truth","h_trig_eff_reco","h_pur_eff_reco"};
  TString xTitle[NTYPE]   = {"#it{p}_{T}^{MC} (GeV/#it{c})",
                             "#it{p}_{T}^{reco,corr} (GeV/#it{c})",
                             "#it{p}_{T}^{reco,corr} (GeV/#it{c})"};

  // Loop over non-ref files
  for (size_t ifi=0; ifi<files.size(); ++ifi) {
    TString fName = files[ifi];
    if (fName == refFile) continue;

    cout << "File: " << fName << endl;

    TFile* f = TFile::Open(fName, "READ");
    if (!f || f->IsZombie()) { cout << "  skip (cannot open)\n"; continue; }

    // loop tag directories from this file (robust if ref missing some tags)
    TIter kIt(f->GetListOfKeys());
    TKey* key = 0;

    while ((key = (TKey*)kIt())) {
      if (std::string(key->GetClassName()) != "TDirectoryFile") continue;
      TString tag = key->GetName();
      if (!tag.Contains("ptlead")) continue;

      for (int t=0; t<NTYPE; ++t) {
        TH1D* hA = GetH(f,    tag, hname[t]);
        TH1D* hB = GetH(fRef, tag, hname[t]);
        if (!hA || !hB) { if(hA)delete hA; if(hB)delete hB; continue; }

        TString baseA = BaseNoExt(fName);
        TString baseB = BaseNoExt(refFile);

        TH1D* hR = MakeRatio(hA, hB, Form("h_ratio_%s__%s", typeKey[t].Data(), tag.Data()));
        delete hA;
        delete hB;

        if (!hR) continue;

        TString pdf = Form("%s/%s/%s/ratio_%s__%s__over__%s.pdf",
                           outDir.Data(), typeKey[t].Data(), tag.Data(),
                           typeKey[t].Data(), baseA.Data(), baseB.Data());

        TString head1 = Form("%s ratio: %s / %s", typeKey[t].Data(), baseA.Data(), baseB.Data());
        TString head2 = tag;

        DrawSave(hR, pdf, xTitle[t], head1, head2);

        delete hR;
      }
    }

    f->Close();
    delete f;
  }

  fRef->Close();
  delete fRef;

  cout << "Done. Output PDFs in: " << outDir << endl;
}