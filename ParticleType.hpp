#ifndef PARTICLE_TYPE_HPP
#define PARTICLE_TYPE_HPP

#include <iostream>

class ParticleType {
public:
  ParticleType(char *name, double mass, int charge);
  char *getName() const;
  double getMass() const;
  int getCharge() const;
  virtual double getWidth() const;
  virtual void Print() const;

private:
  char *const fName;
  double const fMass;
  int const fCharge;
};

#endif
