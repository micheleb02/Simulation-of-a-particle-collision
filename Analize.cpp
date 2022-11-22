#include "Particle.hpp"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"
#include <fstream>

void analize() {

  // Definition of a file containing the histograms and arrays to manage them
  TFile *data = new TFile("Data.root");
  TH1 *HTOT[13];
  TString s[13] = {"h1", "h2", "h3",  "h4",  "h5",  "h6", "h7",
                   "h8", "h9", "h10", "h11", "h7c", "h8c"};
  TString names[7] = {"pion+",   "pion-",   "kaon+", "kaon-",
                      "proton+", "proton-", "kaon*"};

  // Filling the arrays and printing the number of entries
  for (int i = 0; i < 13; ++i) {
    HTOT[i] = (TH1 *)data->Get(s[i]);
  }

  // Defining a fit for both the axis projections of the angles distributions
  // (uniform)
  TF1 *f1 = new TF1("f1", "[0]", 0, 2 * TMath::Pi());
  TF1 *f2 = new TF1("f2", "[0]", 0, TMath::Pi());
  TH1D *ProjX = (((TH2 *)data->Get(s[1]))->ProjectionX("ProjX", 0, 100));
  TH1D *ProjY = (((TH2 *)data->Get(s[1]))->ProjectionY("ProjY", 0, 100));
  ProjX->Fit("f1", "Q0");
  ProjY->Fit("f2", "Q0");

  // Defining a fit for the impulse distribution (exponential)
  TF1 *f3 = new TF1("f3", "[0]*e^(-[1]*x)", 0, 7);
  HTOT[2]->Fit("f3", "Q0");

  // Defining the histograms difference and fitting them with a gaussian
  TH1F *h9c = new TH1F(*(TH1F *)data->Get(s[8]));
  TH1F *h7c = new TH1F(*(TH1F *)data->Get(s[11]));
  h7c->SetTitle("Difference between invariant mass (Depending on charge)");
  TH1F *h11c = new TH1F(*(TH1F *)data->Get(s[10]));

  h7c->Add(HTOT[12], -1);
  h9c->Add(HTOT[9], -1);

  TF1 *f4 = new TF1("f4", "gaus", 0.75, 1.05);
  h7c->Fit("f4", "Q0");

  TF1 *f5 = new TF1("f5", "gaus", 0, 7);
  h9c->Fit("f5", "Q0");

  TF1 *f6 = new TF1("f6", "gaus", 0, 2);
  h11c->Fit("f6", "Q0");

  // Creating a text file containing all results obtained
  std::ofstream fw("./Data result.txt", std::ofstream::out);

  if (fw.is_open()) {

    fw << "##NUMBER OF ENTRIES FOR EACH HISTOGRAM##\n";

    for (int i = 0; i < 11; ++i) {
      fw << HTOT[i]->GetTitle() << " has " << HTOT[i]->GetEntries()
         << " entries" << '\n';
    }

    fw << "\n##PERCENTAGE OF GENERATED PARTICLE##\n";
    for (int i = 0; i < 7; ++i) {
      fw << names[i] << " = ("
         << ((HTOT[0]->GetBinContent(i + 1)) / (1E7)) * 100 << " +/- "
         << (HTOT[0]->GetBinError(i + 1) / (1E7)) * 100 << ") %\n";
    }

    fw << "\n##ANGLES DISTRIBUTION DATA##\n";
    fw << "Fit parameter (phi) = " << f1->GetParameter(0) << " +/- "
       << f1->GetParError(0)
       << "\nFit parameter (theta) = " << f2->GetParameter(0) << " +/- "
       << f2->GetParError(0) << "\nChi squared (phi)= " << f1->GetChisquare()
       << "\nDOF (phi)= " << f1->GetNDF()
       << "\nNormalized Chi (phi) = " << (f1->GetChisquare()) / (f1->GetNDF())
       << "\nChi squared (theta) = " << f2->GetChisquare()
       << "\nDOF (theta) = " << f2->GetNDF()
       << "\nNormalized Chi (theta) = " << (f2->GetChisquare()) / (f2->GetNDF())
       << '\n'
       << '\n';

    fw << "##EXPONENTIAL FIT OF THE IMPULSE DISTRIBUTION##\n";
    fw << "Parameter A = " << f3->GetParameter(0) << " +/- "
       << f3->GetParError(0) << " GeV\nParameter B = " << f3->GetParameter(1)
       << " +/- " << f3->GetParError(1)
       << "\nChi squared = " << f3->GetChisquare() << "\nDOF = " << f3->GetNDF()
       << "\nNormalized Chi = " << ((f3->GetChisquare()) / (f3->GetNDF()))
       << "\nFit probability = " << f3->GetProb() << '\n'
       << '\n';

    fw << "##GAUSSIAN FIT OF THE DIFFERENCE BEWTEEN INVARIANT MASS "
          "DISTRIBUTION "
          "(DEPENDING ON: CHARGE)##\n";
    fw << "Mean (K* mass) = " << f4->GetParameter(1) << " +/- "
       << f4->GetParError(1)
       << " GeV/c^2\nStd. Dev. (K* width) = " << f4->GetParameter(2) << " +/- "
       << f4->GetParError(2) << "GeV/c^2\nAmplitude = " << f4->GetParameter(0)
       << " +/- " << f4->GetParError(0) << " GeV/c^2\nNormalized Chi = "
       << ((f4->GetChisquare()) / (f4->GetNDF()))
       << "\nFit probability = " << f4->GetProb() << '\n'
       << '\n';

    fw << "##GAUSSIAN FIT OF THE DIFFERENCE BEWTEEN INVARIANT MASS "
          "DISTRIBUTION "
          "(DEPENDING ON: TYPE OF PARTICLE)##\n";
    fw << "Mean (K* mass) = " << f5->GetParameter(1) << " +/- "
       << f5->GetParError(1)
       << " GeV/c^2\nStd. Dev. (K* width) = " << f5->GetParameter(2) << " +/- "
       << f5->GetParError(2) << " GeV/c^2\nAmplitude = " << f5->GetParameter(0)
       << " +/- " << f5->GetParError(0) << " GeV/c^2\nNormalized Chi = "
       << ((f5->GetChisquare()) / (f5->GetNDF()))
       << "\nFit probability = " << f5->GetProb() << '\n'
       << '\n';

    fw << "##GAUSSIAN FIT OF THE INVARIANT MASS DISTRIBUTION BETWEEN DECAY "
          "PRODUCTS##\n";
    fw << "Mean = " << f6->GetParameter(1) << " +/- " << f6->GetParError(1)
       << " GeV/c^2\nStd. Dev. = " << f6->GetParameter(2) << " +/- "
       << f6->GetParError(2) << "GeV/c^2\nAmplitude = " << f6->GetParameter(0)
       << " +/- " << f6->GetParError(0) << " GeV/c^2\nNormalized Chi = "
       << ((f6->GetChisquare()) / (f6->GetNDF()))
       << "\nFit probability = " << f6->GetProb() << '\n'
       << '\n';

    fw.close();
  } else {
    std::cout << "Unable tu open Data result.txt" << '\n';
  }

  // Style and canvas for the first pdf + drawing

  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);

  gStyle->SetStatFontSize(0.035);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.08);

  TCanvas *c = new TCanvas("c", "First batch", 200, 10, 1200, 700);
  c->Print("First batch.pdf[");
  c->Divide(2, 2);

  c->cd(1);
  HTOT[0]->GetXaxis()->SetTitle("Type of particle");
  HTOT[0]->GetYaxis()->SetTitle("Number of entries");
  HTOT[0]->GetYaxis()->SetTitleOffset(1.1);
  HTOT[0]->SetTitleSize(0.05, "x");
  HTOT[0]->SetTitleSize(0.05, "y");
  HTOT[0]->SetLineColor(kBlue);
  HTOT[0]->SetFillColor(kAtlantic);
  HTOT[0]->Draw();
  c->cd(2);
  HTOT[2]->GetXaxis()->SetTitle("Impulse");
  HTOT[2]->GetYaxis()->SetTitle("Number of entries");
  HTOT[2]->GetYaxis()->SetTitleOffset(1.1);
  HTOT[2]->SetTitleSize(0.05, "x");
  HTOT[2]->SetTitleSize(0.05, "y");
  HTOT[2]->SetLineColor(kBlue);
  HTOT[2]->SetFillColor(kAtlantic);
  HTOT[2]->Draw();
  f3->Draw("SAME");
  c->cd(3);
  ProjX->SetTitle("Phi angle distribution");
  ProjX->GetXaxis()->SetTitle("Phi");
  ProjX->GetYaxis()->SetTitle("Number of entries");
  ProjX->GetYaxis()->SetTitleOffset(1.1);
  ProjX->SetTitleSize(0.05, "x");
  ProjX->SetTitleSize(0.05, "y");
  ProjX->SetLineColor(kBlue);
  ProjX->SetFillColor(kAtlantic);
  ProjX->Draw();
  f1->Draw("SAME");
  c->cd(4);
  ProjY->SetTitle("Theta angle distribution");
  ProjY->GetXaxis()->SetTitle("Theta");
  ProjY->GetYaxis()->SetTitle("Number of entries");
  ProjY->GetYaxis()->SetTitleOffset(1.1);
  ProjY->SetTitleSize(0.05, "x");
  ProjY->SetTitleSize(0.05, "y");
  ProjY->SetLineColor(kBlue);
  ProjY->SetFillColor(kAtlantic);
  ProjY->Draw();
  f2->Draw("SAME");

  c->Print("First batch.pdf");
  c->Print("First batch.pdf]");

  // Style for the second pdf + drawing

  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.9);

  TCanvas *c2 =
      new TCanvas("c2", "Invariant mass histograms", 200, 10, 400, 800);
  c2->Print("Invariant mass histograms.pdf[");
  c2->Divide(1, 3);

  // Histograms cosmetics
  h11c->GetXaxis()->SetRangeUser(0.6, 1.2);
  h11c->GetYaxis()->SetRangeUser(0., 2500.);
  h7c->GetYaxis()->SetRangeUser(-1200., 6000.);
  h9c->GetXaxis()->SetRangeUser(0.4, 1.4);

  TF1 *FIT[3] = {f4, f5, f6};
  TH1F *HFIT[3] = {h7c, h9c, h11c};

  for (int i = 0; i < 3; ++i) {
    c2->cd(i + 1);
    HFIT[i]->GetXaxis()->SetTitle("Invariant mass (Gev/c^{2})");
    HFIT[i]->SetTitleSize(0.03, "x");
    HFIT[i]->GetYaxis()->SetTitle("Number of entries");
    HFIT[i]->GetYaxis()->SetTitleOffset(1.6);
    HFIT[i]->SetTitleSize(0.03, "y");
    HFIT[i]->SetLabelSize(0.03, "x");
    HFIT[i]->SetLabelSize(0.03, "y");
    HFIT[i]->SetLineColor(kBlue);
    HFIT[i]->SetFillColor(kAtlantic);
    HFIT[i]->Sumw2(kFALSE);
    HFIT[i]->Draw();
    FIT[i]->SetLineColor(kOrange + 7);
    FIT[i]->Draw("SAME");
  }

  c2->Print("Invariant mass histograms.pdf");
  c2->Print("Invariant mass histograms.pdf]");
}