#include "Particle.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>

Particle::Particle(char *name, double Px, double Py, double Pz)
    : fPx{Px}, fPy{Py}, fPz{Pz} {
  fIndex = FindParticle(name);
}

int Particle::fNParticleType = 0;

std::vector<ParticleType *> Particle::fParticleType = {};

int Particle::FindParticle(char *name) {
  int index;
  auto i = std::find_if(fParticleType.begin(), fParticleType.end(),
                        [name](auto temp) {
                          bool check = ((temp->getName()) == name);
                          return check;
                        });
  if (i != fParticleType.end()) {
    index = std::distance(fParticleType.begin(), i);
  } else {
    index = -1;
  }
  return index;
}

int Particle::getIndex() const { return fIndex; }

void Particle::AddParticleType(char *name, double mass, int charge,
                               double width) {

  int size = fParticleType.size();

  ParticleType *newPart;
  if (width == 0.) {
    newPart = new ParticleType(name, mass, charge);
  } else {
    newPart = new ResonanceType(name, mass, charge, width);
  }
  if (size != 0) {
    int index = FindParticle(name);
    if (index == -1) {
      fParticleType.push_back(newPart);
    } else {
      std::cout << "The particle is already there" << '\n';
    }
  } else {
    fParticleType.push_back(newPart);
  }
  ++fNParticleType;
}

void Particle::setter(int index) { fIndex = index; }

void Particle::setter(char *name) { fIndex = FindParticle(name); }

void Particle::PrintArray() {
  for (auto i : fParticleType) {
    i->Print();
  }
}

void Particle::PrintData() const {
  std::cout << "Index= " << fIndex << " ;"
            << "Name= " << *(fParticleType[fIndex]->getName()) << " ;"
            << "(Px, Py, Pz)= "
            << "(" << fPx << "," << fPy << "," << fPz << ")" << '\n';
}

double Particle::getPx() const { return fPx; }

double Particle::getPy() const { return fPy; }

double Particle::getPz() const { return fPz; }

double Particle::getMass() const { return fParticleType[fIndex]->getMass(); }

double Particle::getEnergy() const {
  return std::sqrt(getMass() * getMass() + fPx * fPx + fPy * fPy + fPz * fPz);
}

int Particle::getCharge() const {return fParticleType[getIndex()]->getCharge();}

double Particle::InvMass(Particle &p) const {
  double sumE2 = (getEnergy() + p.getEnergy()) * (getEnergy() + p.getEnergy());
  double sumP2 = (fPx + p.fPx) * (fPx + p.fPx) + (fPy + p.fPy) * (fPy + p.fPy) +
                 (fPz + p.fPz) * (fPz + p.fPz);
  double result = std::sqrt(sumE2 - sumP2);
  return result;
}

void Particle::setP(double px, double py, double pz) {
  fPx = px;
  fPy = py;
  fPz = pz;
}

int Particle::Decay2body(Particle &dau1, Particle &dau2) const {
  if (getMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if (fIndex > -1) { // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;

    massMot += fParticleType[fIndex]->getWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf("Decayment cannot be preformed because mass is too low in this "
           "channel\n");
    return 2;
  }

  double pout =
      sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.setP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi),
            pout * cos(theta));
  dau2.setP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi),
            -pout * cos(theta));

  double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

  double bx = fPx / energy;
  double by = fPy / energy;
  double bz = fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}
void Particle::Boost(double bx, double by, double bz) {

  double energy = getEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / sqrt(1.0 - b2);
  double bp = bx * fPx + by * fPy + bz * fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fPx += gamma2 * bp * bx + gamma * bx * energy;
  fPy += gamma2 * bp * by + gamma * by * energy;
  fPz += gamma2 * bp * bz + gamma * bz * energy;
}