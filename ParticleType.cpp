#include "ParticleType.hpp"

ParticleType::ParticleType(char *name, double mass, int charge)
    : fName{name}, fMass{mass}, fCharge{charge} {}

char *ParticleType::getName() const { return fName; }

double ParticleType::getMass() const { return fMass; }

int ParticleType::getCharge() const { return fCharge; }

double ParticleType::getWidth() const {return 0;}

void ParticleType::Print() const {
  std::cout << "Name= " << *fName << " ,"
            << "Mass= " << fMass << " ,"
            << "Charge= " << fCharge << '\n';
}
