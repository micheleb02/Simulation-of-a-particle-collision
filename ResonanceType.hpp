#ifndef RESONANCE_TYPE_HPP
#define RESONANCE_TYPE_HPP

#include "ParticleType.hpp"

class ResonanceType : public ParticleType {
public:
  ResonanceType(char *name, double mass, int charge, double width);
  double getWidth() const override;
  void Print() const override;

private:
  double const fWidth;
};

#endif
