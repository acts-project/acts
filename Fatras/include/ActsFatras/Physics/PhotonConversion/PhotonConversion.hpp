// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Definitions/Units.hpp"
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
  using Scalar = ActsFatras::Particle::Scalar;
  /// Lower energy threshold for produced children
  Scalar childMomentumMin = 50. * Acts::UnitConstants::MeV;
  /// Scaling factor of children energy
  Scalar childEnergyScaleFactor = 2.;
  /// Scaling factor for photon conversion probability
  Scalar conversionProbScaleFactor = 0.98;

  /// @brief Method for evaluating the distance after which the photon
  /// conversion will occur
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] momentum The momentum of the particle
  ///
  /// @param The distance in X_0
  template <typename generator_t>
  std::pair<Scalar, Scalar> generatePathLimits(generator_t& generator, const Particle& particle) const;

  /// @brief This method evaluates the final state due to the photon conversion
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in, out] particle The interacting photon
  /// @param [out] generated List of generated particles
  ///
  /// @return True if the conversion occured, else false
  template <typename generator_t>
  bool run(generator_t& generator,
                                   Particle& particle, std::vector<Particle>& generated) const;

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
      Scalar childEnergy,
      const Particle::Vector3& childDirection) const;

  /// @brief This method evaluates the energy of a child particle
  ///
  /// @tparam generator_t Type of the random number generator
  /// @param [in, out] generator The random number generator
  /// @param [in] gammaMom The momentum of the photon
  ///
  /// @return The energy of the child particle
  template <typename generator_t>
  Scalar childEnergyFraction(generator_t& generator,
                                       Scalar gammaMom) const;

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
  Scalar screenFunction1(Scalar delta) const;
  Scalar screenFunction2(Scalar delta) const;

  static constexpr Scalar m_Z = 13.;  // Aluminium

  /// Fractions
  static constexpr Scalar m_oneOverThree = 1. / 3.;
  static constexpr Scalar m_nineOverSeven = 9. / 7.;
};

inline Particle::Scalar PhotonConversion::screenFunction1(Scalar delta) const {
  // Compute the value of the screening function 3*PHI1(delta) - PHI2(delta)
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958)
                       : 42.184 - delta * (7.444 - 1.623 * delta);
}

inline Particle::Scalar PhotonConversion::screenFunction2(Scalar delta) const {
  // Compute the value of the screening function 1.5*PHI1(delta)
  // +0.5*PHI2(delta)
  return (delta > 1.4) ? 42.038 - 8.29 * log(delta + 0.958)
                       : 41.326 - delta * (5.848 - 0.902 * delta);
}

template <typename generator_t>
std::pair<Particle::Scalar, Particle::Scalar> PhotonConversion::generatePathLimits(generator_t& generator,
                                                  const Particle& particle) const {
  /// This method is based upon the Athena class PhotonConversionTool

  // Fast exit if not a photon or the energy is too low
  if (particle.pdg() != 22 ||
      particle.absoluteMomentum() < 2. * childMomentumMin)
    return std::make_pair(std::numeric_limits<Scalar>::infinity(), std::numeric_limits<Scalar>::infinity());

  // Use for the moment only Al data - Yung Tsai - Rev.Mod.Particle Physics Vol.
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
  constexpr Scalar p0 = -7.01612e-03;
  constexpr Scalar p1 = 7.69040e-02;
  constexpr Scalar p2 = -6.07682e-01;

  // Calculate xi
  const Scalar xi = p0 + p1 * pow(particle.absoluteMomentum(), p2);

  std::uniform_real_distribution<Scalar> uniformDistribution{0., 1.};
  // This is a transformation of eq. 3.75
  return std::make_pair(-m_nineOverSeven *
          log(conversionProbScaleFactor * (1 - uniformDistribution(generator))) /
          (1. - xi), std::numeric_limits<Scalar>::infinity());
}

template <typename generator_t>
Particle::Scalar PhotonConversion::childEnergyFraction(
    generator_t& generator, Scalar gammaMom) const {
  /// This method is based upon the Geant4 class G4PairProductionRelModel

  /// @note This method is from the Geant4 class G4Element
  //
  //  Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)
  constexpr Scalar k1 = 0.0083;
  constexpr Scalar k2 = 0.20206;
  constexpr Scalar k3 = 0.0020;  // This term is missing in Athena
  constexpr Scalar k4 = 0.0369;
  constexpr Scalar alphaEM = 1. / 137.;
  constexpr Scalar az2 = (alphaEM * m_Z) * (alphaEM * m_Z);
  constexpr Scalar az4 = az2 * az2;
  constexpr Scalar coulombFactor =
      (k1 * az4 + k2 + 1. / (1. + az2)) * az2 - (k3 * az4 + k4) * az4;

  constexpr Scalar logZ13 = log(m_Z) * m_oneOverThree;
  constexpr Scalar FZ = 8. * (logZ13 + coulombFactor);
  const Scalar deltaMax = exp((42.038 - FZ) * 0.1206) - 0.958;

  constexpr Scalar deltaPreFactor = 136. / pow(m_Z, m_oneOverThree);
  const Scalar eps0 = findMass(Acts::eElectron) / gammaMom;
  const Scalar deltaFactor = deltaPreFactor * eps0;
  const Scalar deltaMin = 4. * deltaFactor;

  // Compute the limits of eps
  const Scalar epsMin =
      std::max(eps0, 0.5 - 0.5 * std::sqrt(1. - deltaMin / deltaMax));
  const Scalar epsRange = 0.5 - epsMin;

  // Sample the energy rate (eps) of the created electron (or positron)
  const Scalar F10 = screenFunction1(deltaMin) - FZ;
  const Scalar F20 = screenFunction2(deltaMin) - FZ;
  const Scalar NormF1 = F10 * epsRange * epsRange;
  const Scalar NormF2 = 1.5 * F20;

  // We will need 3 uniform random number for each trial of sampling
  Scalar greject = 0.;
  Scalar eps;
  std::uniform_real_distribution<Scalar> rndmEngine;
  do {
    if (NormF1 > rndmEngine(generator) * (NormF1 + NormF2)) {
      eps = 0.5 - epsRange * pow(rndmEngine(generator), m_oneOverThree);
      const Scalar delta = deltaFactor / (eps * (1. - eps));
      greject = (screenFunction1(delta) - FZ) / F10;
    } else {
      eps = epsMin + epsRange * rndmEngine(generator);
      const Scalar delta = deltaFactor / (eps * (1. - eps));
      greject = (screenFunction2(delta) - FZ) / F20;
    }
  } while (greject < rndmEngine(generator));
  //  End of eps sampling
  return eps * childEnergyScaleFactor;
}

template <typename generator_t>
Particle::Vector3 PhotonConversion::childDirection(
    generator_t& generator,
    const Particle::Vector4& gammaMom4) const {
  /// This method is based upon the Athena class PhotonConversionTool

  // Following the Geant4 approximation from L. Urban
  // the azimutal angle
  Scalar theta = findMass(Acts::eElectron) / gammaMom4[Acts::eEnergy];

  std::uniform_real_distribution<Scalar> uniformDistribution{0., 1.};
  const Scalar u =
      -log(uniformDistribution(generator) * uniformDistribution(generator)) *
      1.6;

  theta *= (uniformDistribution(generator) < 0.25)
               ? u
               : u * m_oneOverThree;  // 9./(9.+27) = 0.25

  // More complex but "more true"
  const Particle::Vector3 gammaMomHep =
      gammaMom4.template segment<3>(Acts::eMom0);
  const Particle::Vector3 newDirectionHep(gammaMomHep.normalized());

  // If it runs along the z axis - no good ==> take the x axis
  const Scalar x =
      (newDirectionHep[Acts::eZ] * newDirectionHep[Acts::eZ] > 0.999999)
          ? 1.
          : -newDirectionHep[Acts::eY];
  const Scalar y = newDirectionHep[Acts::eX];

  // Deflector direction
  const Particle::Vector3 deflectorHep(x, y, 0.);
  // Rotate the new direction for scattering
  Eigen::Transform<Scalar, 3, Eigen::Affine> rotTheta;
  rotTheta = Eigen::AngleAxis<Scalar>(theta, deflectorHep);

  // And arbitrarily in psi
  Eigen::Transform<Scalar, 3, Eigen::Affine> rotPsi;
  rotPsi = Eigen::AngleAxis<Scalar>(
      uniformDistribution(generator) * 2. * M_PI, gammaMomHep);

  return rotPsi * rotTheta * newDirectionHep;
}

template <typename generator_t>
std::vector<Particle> PhotonConversion::recordProduct(
    generator_t& generator, const Particle& photon,
    Scalar childEnergy,
    const Particle::Vector3& childDirection) const {
  // Calculate the child momentum
  const Scalar massChild = findMass(Acts::eElectron);
  const Scalar momentum1 =
      sqrt(childEnergy * childEnergy - massChild * massChild);

  // Use energy-momentum conservation for the other child
  const Particle::Vector3 vtmp =
      photon.fourMomentum().template segment<3>(Acts::eMom0) -
      momentum1 * childDirection;
  const Scalar momentum2 = vtmp.norm();

  // Charge sampling
  std::uniform_int_distribution<> uniformDistribution{0, 1};
  Scalar charge1;
  Scalar charge2;
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
                         .setPosition4(photon.fourPosition())
                         .setDirection(childDirection)
                         .setAbsoluteMomentum(momentum1)
                         .setProcess(ProcessType::ePhotonConversion);
    children.push_back(std::move(child));
  }
  if (momentum2 > childMomentumMin) {
    Particle child = Particle(Barcode(), pdg2)
                         .setPosition4(photon.fourPosition())
                         .setDirection(childDirection)
                         .setAbsoluteMomentum(momentum2)
                         .setProcess(ProcessType::ePhotonConversion);
    children.push_back(std::move(child));
  }
  return children;
}

template <typename generator_t>
bool PhotonConversion::run(
    generator_t& generator, Particle& particle, std::vector<Particle>& generated) const {
	// Fast exit if particle is not alive
  if(!particle)
	return true;
	
	// Fast exit if momentum is too low
  const Scalar p = particle.absoluteMomentum();
  if(p < 2. * childMomentumMin)
	return false;

  // Get one child energy
  const Scalar childEnergy = p * childEnergyFraction(generator, p);

  // Now get the deflection
  const Particle::Vector3 childDir = childDirection(generator, particle.fourMomentum());

  // Produce the final state
  const std::vector<Particle> finalState = recordProduct(generator, particle, childEnergy, childDir);
  generated.insert(generated.end(), finalState.begin(), finalState.end());
  
  // Kill the photon
  particle.setAbsoluteMomentum(0.);
  
  return true;
}
}  // namespace ActsFatras