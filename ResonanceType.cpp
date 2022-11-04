#include "ResonanceType.hpp"

ResonanceType::ResonanceType(char *name, double mass, int charge, double width)
    : ParticleType(name, mass, charge), fWidth{width} {}

double ResonanceType::getWidth() const { return fWidth; }

void ResonanceType::Print() const {
  std::cout << "Name= " << *getName() << " ,"
            << "Mass= " << getMass() << " ,"
            << "Charge= " << getCharge() << " ,"
            << "Width= " << fWidth << '\n';
}
