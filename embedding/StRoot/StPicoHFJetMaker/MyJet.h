#pragma once

#ifndef __CINT__
#include "fastjet/PseudoJet.hh"
#endif

#include <vector>

class MyJet {

public:
  double pt;
  double pt_corr;
  double eta;
  double phi;
  double area;
  double rho;
  double pt_lead;
  int n_constituents;
  double neutral_fraction;
  bool trigger_match;

  MyJet()
      : pt(-9), pt_corr(-9), eta(-9), phi(-9), area(-9), rho(-9), pt_lead(-9),
        n_constituents(-9), neutral_fraction(-9), trigger_match(false) {}

  MyJet(double pt, double pt_corr, double eta, double phi, double area,
        double rho, double pt_lead, int n_constituents, double neutral_fraction,
        bool trigger_match)
      : pt(pt), pt_corr(pt_corr), eta(eta), phi(phi), area(area), rho(rho),
        pt_lead(pt_lead), n_constituents(n_constituents),
        neutral_fraction(neutral_fraction), trigger_match(trigger_match) {}

#ifndef __CINT__
  MyJet(fastjet::PseudoJet jet, double rho)
      : rho(rho) { // instant initialization
    pt = jet.perp();
    eta = jet.eta();
    phi = jet.phi();
    area = jet.area();
    pt_corr = pt - area * rho;
    vector<fastjet::PseudoJet> constituents = sorted_by_pt(jet.constituents());
    pt_lead = constituents[0].perp();
    n_constituents = constituents.size();

    double neutral_sum = 0.0;
    trigger_match = false;
    for (const auto &constituent : constituents) {
      int uidx = constituent.user_index();
      if (uidx == 9999) {
        trigger_match = true;
      }
      if (uidx == 0 || uidx == 9999) { // neutral particles
        neutral_sum += constituent.perp();
      }
    }
    neutral_fraction = neutral_sum / pt;
  }
#endif
};

typedef pair<MyJet, MyJet> MatchedJetPair;
