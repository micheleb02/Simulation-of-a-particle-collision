#include "Particle.hpp"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TStyle.h"

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
  TH2F *h2 = new TH2F("h2", "Angle distribution", 100, 0, 2 * TMath::Pi(), 100,
                      0, TMath::Pi());
  TH1F *h3 = new TH1F("h3", "Impulse distribution", 1000, 0, 5);
  TH1F *h4 = new TH1F("h4", "Trasverse impulse distribution", 1000, 0, 5);
  TH1F *h5 = new TH1F("h5", "Energy distribution", 1000, 0, 10);
  TH1F *h6 = new TH1F("h6", "Invariant mass", 1000, 0, 7);
  h6->Sumw2();
  TH1F *h7 = new TH1F("h7", "Invariant mass - discord charge", 1000, 0, 7);
  h7->Sumw2();
  TH1F *h7c = new TH1F("h7c", "h7 copy to fit", 100, 0.75, 1.05);
  h7c->Sumw2();
  TH1F *h8 =
      new TH1F("h8", "Invariant mass - same charge particles", 1000, 0, 7);
  h8->Sumw2();
  TH1F *h8c = new TH1F("h8c", "h8 copy to fit", 100, 0.75, 1.05);
  h8c->Sumw2();
  TH1F *h9 = new TH1F("h9", "Invariant mass - pion+ & kaon- and pion- & kaon+",
                      1000, 0, 7);
  h9->Sumw2();
  TH1F *h10 = new TH1F(
      "h10", "Invariant mass - pion+ & kaon+ and pion- & kaon-", 1000, 0, 7);
  h10->Sumw2();
  TH1F *h11 = new TH1F("h11", "Invariant mass - decay products", 1000, 0, 2);
  h11->Sumw2();

  TH1 *HTOT[13] = {h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h7c, h8c};

  for (int i = 0; i < 1E5; ++i) {
    for (int j = 0; j < 100; ++j) {

      var = gRandom->Rndm();
      phi = gRandom->Rndm() * TMath::Pi() * 2;
      theta = gRandom->Rndm() * TMath::Pi();

      h2->Fill(phi, theta);

      m = gRandom->Exp(1);

      double px = m * sin(theta) * cos(phi);
      double py = m * sin(theta) * sin(phi);
      double pz = m * cos(theta);

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
        if (a.getPCharge() * b.getPCharge() == -1) {
          h7->Fill(invMass);
          h7c->Fill(invMass);
        }
        if (a.getPCharge() * b.getPCharge() == 1) {
          h8->Fill(invMass);
          h8c->Fill(invMass);
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

  TFile *data = new TFile("Data.root", "RECREATE");

  for (int i = 0; i < 13; ++i) {
    HTOT[i]->Write();
  }

  data->Close();
}