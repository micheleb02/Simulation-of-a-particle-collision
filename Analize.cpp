#include "Particle.hpp"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"

void analize() {
  TFile *data = new TFile("Data.root");
  TH1 *HTOT[11];
  TString s[11] = {"h1", "h2", "h3", "h4",  "h5", "h6",
                   "h7", "h8", "h9", "h10", "h11"};
  TString names[7] = {"pion+",   "pion-",   "kaon+", "kaon-",
                      "proton+", "proton-", "kaon*"};

  for (int i = 0; i < 11; ++i) {
    HTOT[i] = (TH1 *)data->Get(s[i]);
    std::cout << HTOT[i]->GetTitle() << " has  " << HTOT[i]->GetEntries()
              << "  entries" << '\n'
              << '\n';
  }

  std::cout << "Percentage of particles:\n";
  for (int i = 0; i < 7; ++i) {
    std::cout << names[i] << " = " << ((HTOT[0]->GetBinContent(i + 1)) / (1E7))
              << '\n';
  }

  TF1 *f1 = new TF1("f1", "[0]", 0, 2 * TMath::Pi());
  TF1 *f2 = new TF1("f2", "[0]", 0, TMath::Pi());
  TH1D *ProjX = (((TH2 *)data->Get(s[1]))->ProjectionX("ProjX", 0, 100));
  TH1D *ProjY = (((TH2 *)data->Get(s[1]))->ProjectionX("ProjY", 0, 100));
  ProjX->Fit("f1", "Q0");
  ProjY->Fit("f2", "Q0");

  std::cout << "\nMean phi = " << f1->GetParameter(0)
            << "\nMean theta = " << f2->GetParameter(0) << '\n'
            << '\n';
  std::cout << "Normalized Chi for phi = "
            << (f1->GetChisquare()) / (f1->GetNDF())
            << "\nNormalized Chi for theta  = "
            << (f2->GetChisquare()) / (f2->GetNDF()) << '\n'
            << '\n';

  TF1 *f3 = new TF1("f3", "[0]*e^([1]*x)", 0, 7);
  HTOT[2]->Fit("f3", "Q0");
  std::cout << "Parameter A =" << f3->GetParameter(0)
            << "\nParameter B = " << f3->GetParameter(1)
            << "\nNormalized Chi = " << ((f3->GetChisquare()) / (f3->GetNDF()))
            << "\nFit probability = " << f3->GetProb() << '\n'
            << '\n';
}