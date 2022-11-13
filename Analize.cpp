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
    std::cout << HTOT[i]->GetTitle() << " has " << HTOT[i]->GetEntries()
              << " entries" << '\n';
  }

  std::cout << "\nPercentage of particles:\n";
  for (int i = 0; i < 7; ++i) {
    std::cout << names[i] << " = " << ((HTOT[0]->GetBinContent(i + 1)) / (1E7))
              << '\n';
  }

  TF1 *f1 = new TF1("f1", "[0]", 0, 2 * TMath::Pi());
  TF1 *f2 = new TF1("f2", "[0]", 0, TMath::Pi());
  TH1D *ProjX = (((TH2 *)data->Get(s[1]))->ProjectionX("ProjX", 0, 100));
  TH1D *ProjY = (((TH2 *)data->Get(s[1]))->ProjectionY("ProjY", 0, 100));
  ProjX->Fit("f1", "Q0");
  ProjY->Fit("f2", "Q0");

  double meanPhi = 0.;
  double meanTheta = 0.;

  for (int i = 1; i <= 100; ++i) {
    double p = i * (2 * TMath::Pi() / 100.);
    double t = i * (TMath::Pi() / 100.);
    meanPhi += p * (f1->GetParameter(0));
    meanTheta += t * (f2->GetParameter(0));
  }

  meanPhi = meanPhi / (ProjX->GetEntries());
  meanTheta = meanTheta / (ProjY->GetEntries());

  std::cout << "\nMean phi = " << meanPhi << "\nMean theta = " << meanTheta;
  std::cout << "\nNormalized Chi for phi = "
            << (f1->GetChisquare()) / (f1->GetNDF())
            << "\nNormalized Chi for theta  = "
            << (f2->GetChisquare()) / (f2->GetNDF()) << '\n'
            << '\n';

  TF1 *f3 = new TF1("f3", "[0]*e^([1]*x)", 0, 7);
  HTOT[2]->Fit("f3", "Q0");
  std::cout << "Parameter A = " << f3->GetParameter(0) << " +/- "
            << f3->GetParError(0) << "\nParameter B = " << f3->GetParameter(1)
            << " +/- " << f3->GetParError(1)
            << "\nNormalized Chi = " << ((f3->GetChisquare()) / (f3->GetNDF()))
            << "\nFit probability = " << f3->GetProb() << '\n'
            << '\n';

  TH1F *h12 = new TH1F(*(TH1F *)data->Get(s[6]));
  TH1F *h13 = new TH1F(*(TH1F *)data->Get(s[8]));
  h12->Add(h12, HTOT[7], 1, -1);
  h13->Add(h13, HTOT[9], 1, -1);

  TF1 *f4 = new TF1("f4", "gaus", 0, 10);
  h12->Fit("f4", "Q0");
  std::cout << "Mean (K* mass) = " << f4->GetParameter(1) << " +/- "
            << f4->GetParError(1)
            << "\nStd. Dev. (K* width) = " << f4->GetParameter(2) << " +/- "
            << f4->GetParError(2)
            << "\nNormalized Chi = " << ((f4->GetChisquare()) / (f4->GetNDF()))
            << "\nFit probability = " << f4->GetProb() << '\n'
            << '\n';

  TF1 *f5 = new TF1("f5", "gaus", 0, 10);
  h13->Fit("f5", "Q0");
  std::cout << "Mean (K* mass) = " << f5->GetParameter(1) << " +/- "
            << f5->GetParError(1)
            << "\nStd. Dev. (K* width) = " << f5->GetParameter(2) << " +/- "
            << f5->GetParError(2)
            << "\nNormalized Chi = " << ((f5->GetChisquare()) / (f5->GetNDF()))
            << "\nFit probability = " << f5->GetProb() << '\n'
            << '\n';
}