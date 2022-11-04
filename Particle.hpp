#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "ResonanceType.hpp"
#include <vector>

class Particle {
public:
  Particle(char *name, double Px = 0., double Py = 0., double Pz = 0.);
  int getIndex() const;
  void static AddParticleType(char *name, double mass, int charge,
                              double width = 0.);
  void setter(int index);
  void setter(char *name);
  static void PrintArray();
  void PrintData() const;
  double getPx() const;
  double getPy() const;
  double getPz() const;
  double getMass() const;
  double getEnergy() const;
  int getPCharge() const;
  double InvMass(Particle &p) const;
  void setP(double px, double py, double pz);
  int Decay2body(Particle &dau1, Particle &dau2) const;

private:
  static int const fMaxNumParticleType = 10;
  static std::vector<ParticleType *> fParticleType;
  static int fNParticleType;
  int fIndex;
  double fPx;
  double fPy;
  double fPz;
  static int FindParticle(char *name);
  void Boost(double bx, double by, double bz);
};

#endif
