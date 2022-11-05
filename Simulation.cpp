#include "Particle.hpp"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TRandom.h"

void simulation() {

  gRandom->SetSeed();

  std::vector<Particle> EventParticle;
  std::vector<Particle> Resonance;
  std::vector<Particle> Decay;

  char *p = new char{'p'};
  char *P = new char{'P'};
  char *k = new char{'k'};
  char *K = new char{'K'};
  char *d = new char{'d'};
  char *D = new char{'D'};
  char *n = new char{'n'};

  Particle::AddParticleType(p, 0.13597, 1);
  Particle::AddParticleType(P, 0.13597, -1);
  Particle::AddParticleType(k, 0.49267, 1);
  Particle::AddParticleType(K, 0.49267, -1);
  Particle::AddParticleType(d, 0.93827, 1);
  Particle::AddParticleType(D, 0.93827, -1);
  Particle::AddParticleType(n, 0.89166, 0, 0.050);

  double var = 0;
  double theta = 0;
  double phi = 0;
  double m = 0;

  TH1F *h1 = new TH1F("h1", "Type distribution", 7, 0, 7);
  TH3F *h2 =
      new TH3F("h2", "Angle distribution", 100, -1, 1, 100, -1, 1, 100, -1, 1);
  TH1F *h3 = new TH1F("h3", "Impulse distribution", 1000, 0, 5);
  TH1F *h4 = new TH1F("h4", "Trasverse impulse distribution", 1000, 0, 5);
  TH1F *h5 = new TH1F("h5", "Energy distribution", 1000, 0, 10);
  TH1F *h6 = new TH1F("h6", "Invariant mass", 1000, 0, 10);
  h6->Sumw2();
  TH1F *h7 =
      new TH1F("h7", "Invariant mass between discord charge", 1000, 0, 10);
  h7->Sumw2();
  TH1F *h8 = new TH1F("h8", "Invariant mass between same charge particles",
                      1000, 0, 10);
  h8->Sumw2();
  TH1F *h9 = new TH1F("h9", "Invariant mass between specific particles (1)",
                      1000, 0, 10);
  h9->Sumw2();
  TH1F *h10 = new TH1F("h10", "Invariant mass between specific particles (2)",
                       1000, 0, 10);
  h10->Sumw2();
  TH1F *h11 =
      new TH1F("h11", "Invariant mass between decayed products", 1000, 0, 10);
  h11->Sumw2();

  TH1 *HTOT[11] = {h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11};

  TCanvas *c1 = new TCanvas("c1", "Type distribution", 200, 10, 600, 400);
  TCanvas *c2 = new TCanvas("c2", "Angle distribution", 200, 10, 600, 400);
  TCanvas *c3 = new TCanvas("c3", "Impulse distribution", 200, 10, 600, 400);
  TCanvas *c4 =
      new TCanvas("c4", "Trasverse impulse distribution", 200, 10, 600, 400);
  TCanvas *c5 = new TCanvas("c5", "Energy distribution", 200, 10, 600, 400);
  TCanvas *c6 = new TCanvas("c6", "Invariant mass", 200, 10, 600, 400);
  TCanvas *c7 = new TCanvas("c7", "Invariant mass between discord charge", 200,
                            10, 600, 400);
  TCanvas *c8 = new TCanvas(
      "c8", "Invariant mass between same charge particles", 200, 10, 600, 400);
  TCanvas *c9 = new TCanvas(
      "c9", "Invariant mass between specific particles (1)", 200, 10, 600, 400);
  TCanvas *c10 =
      new TCanvas("c10", "Invariant mass between specific particles (2)", 200,
                  10, 600, 400);
  TCanvas *c11 = new TCanvas("c11", "Invariant mass between decayed products",
                             200, 10, 600, 400);

  TCanvas *CTOT[11] = {c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11};

  for (int i = 0; i < 1E5; ++i) {
    for (int j = 0; j < 100; ++j) {

      var = gRandom->Rndm();
      theta = gRandom->Rndm() * TMath::Pi() * 2;
      phi = gRandom->Rndm() * TMath::Pi();

      h2->Fill(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));

      m = gRandom->Exp(1);

      double px = m * sin(phi) * cos(theta);
      double py = m * sin(phi) * sin(theta);
      double pz = m * cos(phi);

      h3->Fill(m);
      h4->Fill(sqrt(px * px + py * py));

      if (var < 0.4) {
        Particle a = Particle(p, px, py, pz);
        EventParticle.push_back(a);
      } else if (var >= 0.4 && var < 0.8) {
        Particle a = Particle(P, px, py, pz);
        EventParticle.push_back(a);
      } else if (var >= 0.8 && var < 0.85) {
        Particle a = Particle(k, px, py, pz);
        EventParticle.push_back(a);
      } else if (var >= 0.85 && var < 0.9) {
        Particle a = Particle(K, px, py, pz);
        EventParticle.push_back(a);
      } else if (var >= 0.9 && var < 0.945) {
        Particle a = Particle(d, px, py, pz);
        EventParticle.push_back(a);
      } else if (var >= 0.945 && var < 0.99) {
        Particle a = Particle(D, px, py, pz);
        EventParticle.push_back(D);
      } else if (var >= 0.99 && var < 0.995) {
        Particle a = Particle(n, px, py, pz);
        Particle b = Particle(P);
        Particle c = Particle(k);
        a.Decay2body(b, c);
        Decay.push_back(a);
        Resonance.push_back(b);
        Resonance.push_back(c);
        h11->Fill(b.InvMass(c));
      } else {
        Particle a = Particle(n, px, py, pz);
        Particle b = Particle(p);
        Particle c = Particle(K);
        a.Decay2body(b, c);
        Decay.push_back(a);
        Resonance.push_back(b);
        Resonance.push_back(c);
        h11->Fill(b.InvMass(c));
      }
    }

    for (auto k : EventParticle) {
      h1->Fill(k.getIndex());
      h5->Fill(k.getEnergy());
    }

    for (auto k : Decay) {
      h1->Fill(k.getIndex());
      h5->Fill(k.getEnergy());
    }

    for (auto k : Resonance) {
      EventParticle.push_back(k);
    }

    int size = EventParticle.size();

    for (int i = 0; i < size - 1; ++i) {
      for (int j = i + 1; j < size; ++j) {
        Particle a = EventParticle[i];
        Particle b = EventParticle[j];
        double invMass = a.InvMass(b);
        h6->Fill(invMass);
        if (a.getPCharge() != b.getPCharge()) {
          h7->Fill(invMass);
        }
        if (a.getPCharge() == b.getPCharge()) {
          h8->Fill(invMass);
        }
        if ((a.getIndex() == 0 && b.getIndex() == 3) ||
            (a.getIndex() == 1 && b.getIndex() == 2)) {
          h9->Fill(invMass);
        }
        if ((a.getIndex() == 0 && b.getIndex() == 2) ||
            (a.getIndex() == 1 && b.getIndex() == 3)) {
          h10->Fill(invMass);
        }
      }
    }

    EventParticle.clear();
    Resonance.clear();
    Decay.clear();
  }

  TFile *Histos = new TFile("Particles_histo.root", "RECREATE");
  for (int i = 0; i < 11; ++i) {
    CTOT[i]->cd();
    HTOT[i]->DrawCopy();
    HTOT[i]->Write();
  }

  Histos->Close();
  
}