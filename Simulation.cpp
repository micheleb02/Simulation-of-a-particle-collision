#include "Particle.hpp"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TRandom.h"

void simulation() {

  gRandom->SetSeed();

  std::vector<Particle> EventParticle;
  std::vector<Particle> Resonance;

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
  TH1F *h7 =
      new TH1F("h7", "Invariant mass between discord charge", 1000, 0, 10);
  TH1F *h8 = new TH1F("h8", "Invariant mass between specific particles (1)",
                      1000, 0, 10);
  TH1F *h9 = new TH1F("h9", "Invariant mass between specific particles (2)",
                      1000, 0, 10);
  TH1F *h10 =
      new TH1F("h10", "Invariant mass between decayed products", 1000, 0, 10);

  TCanvas *c1 = new TCanvas("c1", "Type distribution", 200, 10, 600, 400);
  TCanvas *c2 = new TCanvas("c2", "Angle distribution", 200, 10, 600, 400);
  TCanvas *c3 = new TCanvas("c3", "Impulse distribution", 200, 10, 600, 400);
  TCanvas *c4 =
      new TCanvas("c4", "Trasverse impulse distribution", 200, 10, 600, 400);
  TCanvas *c5 = new TCanvas("c5", "Energy distribution", 200, 10, 600, 400);
  TCanvas *c6 = new TCanvas("c6", "Invariant mass", 200, 10, 600, 400);

  for (int i = 0; i < 10000; ++i) {
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
        EventParticle.push_back(a);
        Resonance.push_back(b);
        Resonance.push_back(c);
      } else {
        Particle a = Particle(n, px, py, pz);
        Particle b = Particle(p);
        Particle c = Particle(K);
        a.Decay2body(b, c);
        EventParticle.push_back(a);
        Resonance.push_back(b);
        Resonance.push_back(c);
      }
    }

    for (auto k : EventParticle) {
      h1->Fill(k.getIndex());
      h5->Fill(k.getEnergy());
    }

    for (auto k : Resonance) {
      EventParticle.push_back(k);
    }

    int size = EventParticle.size();
    int nDecay = Resonance.size();

    for (int i = 0; i < size; ++i) {
      for (int j = i + 1; j < size; ++j) {
        Particle a = EventParticle[i];
        Particle b = EventParticle[j];
        double invMass = a.InvMass(b);
        h6->Fill(invMass);
       /* if (a.getPCharge() != b.getPCharge()) {
          h7->Fill(invMass);
          continue;
        }*/
       /* if (a.getPCharge() == b.getPCharge()) {
          h7->Fill(invMass);
          continue;
        }*/
        if ((a.getIndex() == 0 && b.getIndex() == 3) ||
            (a.getIndex() == 1 && b.getIndex() == 2)) {
          h8->Fill(invMass);
          continue;
        }
        if ((a.getIndex() == 0 && b.getIndex() == 2) ||
            (a.getIndex() == 1 && b.getIndex() == 3)) {
          h9->Fill(invMass);
          continue;
        }
        if (i > (size - nDecay)) {
          h10->Fill(invMass);
          continue;
        }
      }
    }

    EventParticle.clear();
    Resonance.clear();
  }

  c1->cd();
  h1->Draw();
  c2->cd();
  h2->Draw();
  c3->cd();
  h3->Draw();
  c4->cd();
  h4->Draw();
  c5->cd();
  h5->Draw();
  c6->cd();
  h10->Draw();
}