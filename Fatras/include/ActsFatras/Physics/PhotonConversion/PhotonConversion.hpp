// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"
#include "ActsFatras/Utilities/ParticleData.hpp"
#include <limits>
#include <random>
#include <vector>
#include <math.h>

namespace ActsFatras {
/// @brief This class handles the photon conversion. It evaluates the distance
/// after which the interaction will occur and the final state due the
/// interaction itself.
class PhotonConversion {
 public:
  /// Lower energy threshold for produced children
  double childMomentumMin = 50. * Acts::UnitConstants::MeV;
  /// Scaling factor of children energy
  double childEnergyScaleFactor = 2.;
  /// Scaling factor for photon conversion probability
  double conversionProbScaleFactor = 0.98;

  /// @brief Method for evaluating the distance after which the photon
  /// conversion will occur
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] momentum The momentum of the particle
  ///
  /// @param The distance in X_0
  template <typename generator_t>
  void pairProduction(generator_t& generator, Particle& particle) const;

  /// @brief This method evaluates the final state due to the photon conversion
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] particle The interacting photon
  ///
  /// @return Vector containing the final state leptons
  template <typename generator_t>
  std::vector<Particle> operator()(generator_t& generator,
                                   const Acts::MaterialSlab& /*slab*/,
                                   const Particle& particle) const;

 private:
  /// @brief This method constructs and returns the child particles
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] photon The interacting photon
  /// @param [in] childEnergy The energy of one child particle
  /// @param [in] childDirection The direction of the child particle
  ///
  /// @return Vector containing the produced leptons
  template <typename generator_t>
  std::vector<Particle> recordProduct(
      generator_t& generator, const Particle& photon,
      Particle::Scalar childEnergy,
      const Particle::Vector3& childDirection) const;

  /// @brief This method evaluates the energy of a child particle
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] gammaMom The momentum of the photon
  ///
  /// @return The energy of the child particle
  template <typename generator_t>
  Particle::Scalar childEnergyFraction(generator_t& generator,
                                       Particle::Scalar gammaMom) const;

  /// @brief This method evaluates the direction of a child particle
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] gammaMom4 The momentum four vector of the photon
  ///
  /// @return The direction vector of the child particle
  template <typename generator_t>
  Particle::Vector3 childDirection(generator_t& generator,
                                   const Particle::Vector4& gammaMom4) const;

  /// Helper methods for momentum evaluation
  /// @note These methods are taken from the Geant4 class
  /// G4PairProductionRelModel
  double screenFunction1(double delta) const;
  double screenFunction2(double delta) const;

  static constexpr float m_Z = 13.;  // Aluminium

  /// Irrational numbers
  static constexpr double m_oneOverThree = 1. / 3.;
  static constexpr double m_nineOverSeven = 9 / 7;

  //~ constexpr double computeCoulombFactor() const;
};

inline double PhotonConversion::screenFunction1(double delta) const {
  // Compute the value of the screening function 3*PHI1(delta) - PHI2(delta)
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958)
                       : 42.184 - delta * (7.444 - 1.623 * delta);
}

inline double PhotonConversion::screenFunction2(double delta) const {
  // Compute the value of the screening function 1.5*PHI1(delta)
  // +0.5*PHI2(delta)
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958)
                       : 41.326 - delta * (5.848 - 0.902 * delta);
}

//~ inline constexpr double
//~ PhotonConversion::computeCoulombFactor() const
//~ {

//~ }
}  // namespace ActsFatras

template <typename generator_t>
void ActsFatras::PhotonConversion::pairProduction(generator_t& generator,
                                                  Particle& particle) const {
  /// This method is based upon the Athena class PhotonConversionTool

  // Fast exit if not a photon or the energy is too low
  if (particle.pdg() != 22 ||
      particle.absMomentum() < 100. * Acts::UnitConstants::MeV)
    return;

  // use for the moment only Al data - Yung Tsai - Rev.Mod.Particle Physics Vol.
  // 46, No.4, October 1974 optainef from a fit given in the momentum range 100
  // 10 6 2 1 0.6 0.4 0.2 0.1 GeV

  //// Quadratic background function
  //  Double_t fitFunction(Double_t *x, Double_t *par) {
  //  return par[0] + par[1]*pow(x[0],par[2]);
  // }

  // EXT PARAMETER                                   STEP         FIRST
  // NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
  //  1  p0          -7.01612e-03   8.43478e-01   1.62766e-04   1.11914e-05
  //  2  p1           7.69040e-02   1.00059e+00   8.90718e-05  -8.41167e-07
  //  3  p2          -6.07682e-01   5.13256e+00   6.07228e-04  -9.44448e-07
  constexpr double p0 = -7.01612e-03;
  constexpr double p1 = 7.69040e-02;
  constexpr double p2 = -6.07682e-01;

  // calculate xi
  const double xi = p0 + p1 * pow(particle.absMomentum(), p2);

  /// Athena
  //~ double attenuation = exp(
  //-7.777e-01*pathCorrection*mprop.thicknessInX0()*(1.-xi) ); ~ return
  //(conversionProbScaleFactor*CLHEP::RandFlat::shoot(m_randomEngine) >
  // attenuation) ? true : false;

  /// Transformation steps
  //~ return (CLHEP::RandFlat::shoot(m_randomEngine) < 1 - attenuation /
  // conversionProbScaleFactor) ? true : false; ~ conversionProbScaleFactor(rnd
  // - 1) < -attenuation ~ conversionProbScaleFactor(1 - rnd) > attenuation ~
  // ln(conversionProbScaleFactor(1 - rnd)) >
  //-7/9*pathCorrection*mprop.thicknessInX0()*(1.-xi)

  std::uniform_real_distribution<double> uniformDistribution{0., 1.};
  // This is a transformation of eq. 3.75
  particle.setMaterialLimits(
      -m_nineOverSeven *
          ln(conversionProbScaleFactor * (1 - uniformDistribution(generator))) /
          (1. - xi),
      particle.pathLimitL0());
}

template <typename generator_t>
ActsFatras::Particle::Scalar ActsFatras::PhotonConversion::childEnergyFraction(
    generator_t& generator, Particle::Scalar gammaMom) const {
  /// This method is based upon the Geant4 class G4PairProductionRelModel

  /// @note This method is from the Geant4 class G4Element
  //
  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)
  constexpr double k1 = 0.0083;
  constexpr double k2 = 0.20206;
  constexpr double k3 = 0.0020;  // This term is missing in Athena
  constexpr double k4 = 0.0369;
  constexpr double alphaEM = 1. / 137.;
  constexpr double az2 = (alphaEM * m_Z) * (alphaEM * m_Z);
  constexpr double az4 = az2 * az2;
  constexpr double coulombFactor =
      (k1 * az4 + k2 + 1. / (1. + az2)) * az2 - (k3 * az4 + k4) * az4;

  constexpr double logZ13 = log(m_Z) * m_oneOverThree;
  constexpr double FZ = 8. * (logZ13 + coulombFactor);
  const double deltaMax = exp((42.038 - FZ) * 0.1206) - 0.958;

  constexpr double deltaPreFactor = 136. / pow(m_Z, m_oneOverThree);
  const double eps0 = findMass(Acts::eElectron) / gammaMom;
  const double deltaFactor = deltaPreFactor * eps0;
  const double deltaMin = 4. * deltaFactor;

  // compute the limits of eps
  const double epsMin =
      std::max(eps0, 0.5 - 0.5 * std::sqrt(1. - deltaMin / deltaMax));
  const double epsRange = 0.5 - epsMin;

  // sample the energy rate (eps) of the created electron (or positron)
  const double F10 = screenFunction1(deltaMin) - FZ;
  const double F20 = screenFunction2(deltaMin) - FZ;
  const double NormF1 = F10 * epsRange * epsRange;
  const double NormF2 = 1.5 * F20;

  // we will need 3 uniform random number for each trial of sampling
  double greject = 0.;
  double eps;
  std::uniform_real_distribution<double> rndmEngine;
  do {
    if (NormF1 > rndmEngine(generator) * (NormF1 + NormF2)) {
      eps = 0.5 - epsRange * pow(rndmEngine(generator), m_oneOverThree);
      const double delta = deltaFactor / (eps * (1. - eps));
      greject = (screenFunction1(delta) - FZ) / F10;
    } else {
      eps = epsMin + epsRange * rndmEngine(generator);
      const double delta = deltaFactor / (eps * (1. - eps));
      greject = (screenFunction2(delta) - FZ) / F20;
    }
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while (greject < rndmEngine(generator));
  //  end of eps sampling
  return eps * childEnergyScaleFactor;
}

template <typename generator_t>
ActsFatras::Particle::Vector3 ActsFatras::PhotonConversion::childDirection(
    generator_t& generator,
    const ActsFatras::Particle::Vector4& gammaMom4) const {
  /// This method is based upon the Athena class PhotonConversionTool

  // Following the Geant4 approximation from L. Urban
  // the azimutal angle

  // the start of the equation
  Particle::Scalar theta = findMass(Acts::eElectron) / gammaMom4[Acts::eEnergy];

  std::uniform_real_distribution<Particle::Scalar> uniformDistribution{0., 1.};
  const Particle::Scalar u =
      -log(uniformDistribution(generator) * uniformDistribution(generator)) *
      1.6;

  theta *= (uniformDistribution(generator) < 0.25)
               ? u
               : u * m_oneOverThree;  // 9./(9.+27) = 0.25

  // more complex but "more true"
  const Particle::Vector3 gammaMomHep =
      gammaMom4.template segment<3>(Acts::eMom0);
  const Particle::Vector3 newDirectionHep(gammaMomHep.normalized());

  // if it runs along the z axis - no good ==> take the x axis
  const Particle::Scalar x =
      (newDirectionHep[Acts::eZ] * newDirectionHep[Acts::eZ] > 0.999999)
          ? 1.
          : -newDirectionHep[Acts::eY];
  const Particle::Scalar y = newDirectionHep[Acts::eX];

  // deflector direction
  const Particle::Vector3 deflectorHep(x, y, 0.);
  // rotate the new direction for scattering
  Eigen::Transform<Particle::Scalar, 3, Eigen::Affine> rotTheta;
  rotTheta = Eigen::AngleAxis<Particle::Scalar>(theta, deflectorHep);

  // and arbitrarily in psi
  Eigen::Transform<Particle::Scalar, 3, Eigen::Affine> rotPsi;
  rotPsi = Eigen::AngleAxis<Particle::Scalar>(
      uniformDistribution(generator) * 2. * M_PI, gammaMomHep);

  return rotPsi * rotTheta * newDirectionHep;
}

template <typename generator_t>
std::vector<ActsFatras::Particle> ActsFatras::PhotonConversion::recordProduct(
    generator_t& generator, const ActsFatras::Particle& photon,
    ActsFatras::Particle::Scalar childEnergy,
    const ActsFatras::Particle::Vector3& childDirection) const {
  // Calculate the child momentum
  const float massChild = findMass(Acts::eElectron);
  const Particle::Scalar momentum1 =
      sqrt(childEnergy * childEnergy - massChild * massChild);

  // Use energy-momentum conservation for the other child
  const Particle::Vector3 vtmp =
      photon.momentum4().template segment<3>(Acts::eMom0) -
      momentum1 * childDirection;
  const Particle::Scalar momentum2 = vtmp.norm();

  // Charge sampling
  std::uniform_int_distribution<> uniformDistribution{0, 1};
  Particle::Scalar charge1;
  Particle::Scalar charge2;
  if (uniformDistribution(generator) == 1) {
    charge1 = -1.;
    charge2 = 1.;
  } else {
    charge1 = 1.;
    charge2 = -1.;
  }

  // Assign the PDG ID
  const Acts::PdgParticle pdg1 =
      static_cast<Acts::PdgParticle>(std::copysign(Acts::eElectron, charge1));
  const Acts::PdgParticle pdg2 =
      static_cast<Acts::PdgParticle>(std::copysign(Acts::eElectron, charge2));

  // Build the particles and store them
  std::vector<Particle> children;
  children.reserve(2);

  if (momentum1 > childMomentumMin) {
    Particle child = Particle(Barcode(), pdg1)
                         .setPosition4(photon.position4())
                         .setDirection(childDirection)
                         .setAbsMomentum(momentum1)
                         .setProcess(ProcessType::ePhotonConversion);
    children.push_back(std::move(child));
  }
  if (momentum2 > childMomentumMin) {
    Particle child = Particle(Barcode(), pdg2)
                         .setPosition4(photon.position4())
                         .setDirection(childDirection)
                         .setAbsMomentum(momentum2)
                         .setProcess(ProcessType::ePhotonConversion);
    children.push_back(std::move(child));
  }
  return children;
}

template <typename generator_t>
std::vector<ActsFatras::Particle> ActsFatras::PhotonConversion::operator()(
    generator_t& generator, const Acts::MaterialSlab& /*slab*/,
    const ActsFatras::Particle& particle) const {
  const double p = particle.absMomentum();

  // Get one child energy
  const Particle::Scalar childEnergy = p * childEnergyFraction(generator, p);

  // Now get the deflection
  Particle::Vector3 childDir = childDirection(generator, particle.momentum4());

  // TODO: The parent must die, maybe here?

  // Produce the final state
  return recordProduct(generator, particle, childEnergy, childDir);
}