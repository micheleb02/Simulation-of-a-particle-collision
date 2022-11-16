#include "Particle.hpp"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"

void analize() {

  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  TFile *data = new TFile("Data.root");
  TH1 *HTOT[13];
  TString s[13] = {"h1", "h2", "h3",  "h4",  "h5",  "h6", "h7",
                   "h8", "h9", "h10", "h11", "h7c", "h8c"};
  TString names[7] = {"pion+",   "pion-",   "kaon+", "kaon-",
                      "proton+", "proton-", "kaon*"};

  for (int i = 0; i < 13; ++i) {
    HTOT[i] = (TH1 *)data->Get(s[i]);
    if (i < 11) {
      std::cout << HTOT[i]->GetTitle() << " has " << HTOT[i]->GetEntries()
                << " entries" << '\n';
    }
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

  std::cout << "\nMean phi = " << ProjX->GetMean() << " +/- "
            << ProjX->GetMeanError() << "\nMean theta = " << ProjY->GetMean()
            << " +/- " << ProjY->GetMeanError();
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

  TH1F *h9c = new TH1F(*(TH1F *)data->Get(s[8]));
  HTOT[11]->Add(HTOT[12], -1);
  h9c->Add(HTOT[9], -1);

  TF1 *f4 = new TF1("f4", "gaus", 0.75, 1.05);
  HTOT[11]->Fit("f4", "Q0");
  std::cout << "Mean (K* mass) = " << f4->GetParameter(1) << " +/- "
            << f4->GetParError(1)
            << "\nStd. Dev. (K* width) = " << f4->GetParameter(2) << " +/- "
            << f4->GetParError(2)
            << "\nNormalized Chi = " << ((f4->GetChisquare()) / (f4->GetNDF()))
            << "\nFit probability = " << f4->GetProb() << '\n'
            << '\n';

  TF1 *f5 = new TF1("f5", "gaus", 0, 10);
  h9c->Fit("f5", "Q0");
  std::cout << "Mean (K* mass) = " << f5->GetParameter(1) << " +/- "
            << f5->GetParError(1)
            << "\nStd. Dev. (K* width) = " << f5->GetParameter(2) << " +/- "
            << f5->GetParError(2)
            << "\nNormalized Chi = " << ((f5->GetChisquare()) / (f5->GetNDF()))
            << "\nFit probability = " << f5->GetProb() << '\n'
            << '\n';

  TCanvas *c = new TCanvas();
  c->Print("Invariant mass histograms.pdf[");
  c->Divide(2, 3);

  for (int i = 1; i <= 6; ++i) {
    c->cd(i);
    HTOT[i + 4]->GetXaxis()->SetTitle("Invariant mass (Gev/c^{2})");
    HTOT[i + 4]->GetYaxis()->SetTitle("Number of entries");
    HTOT[i + 4]->Sumw2(kFALSE);
    HTOT[i + 4]->Draw();
  }

  c->Print("Invariant mass histograms.pdf");
  c->Print("Invariant mass histograms.pdf]");
}