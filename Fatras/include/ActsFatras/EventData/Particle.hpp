// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <iosfwd>
#include <limits>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

namespace ActsFatras {

/// Simulation particle information and kinematic state.
class Particle {
 public:
  using Scalar = double;
  using Vector3 = Acts::ActsVector<Scalar, 3>;
  using Vector4 = Acts::ActsVector<Scalar, 4>;

  /// Construct a default particle with invalid identity.
  Particle() = default;
  /// Construct a particle at rest with explicit mass and charge.
  ///
  /// @param particleId Particle identifier within an event
  /// @param pdg PDG id
  /// @param charge Particle charge in native units
  /// @param mass Particle mass in native units
  ///
  /// @warning It is the users responsibility that charge and mass match
  ///          the PDG particle number.
  Particle(Barcode particleId, Acts::PdgParticle pdg, Scalar charge,
           Scalar mass)
      : m_particleId(particleId), m_pdg(pdg), m_charge(charge), m_mass(mass) {}
  /// Construct a particle at rest from a PDG particle number.
  ///
  /// @param particleId Particle identifier within an event
  /// @param pdg PDG particle number
  ///
  /// Charge and mass are retrieved from the particle data table.
  Particle(Barcode particleId, Acts::PdgParticle pdg);
  Particle(const Particle &) = default;
  Particle(Particle &&) = default;
  Particle &operator=(const Particle &) = default;
  Particle &operator=(Particle &&) = default;

  /// Construct a new particle with a new identifier but same kinematics.
  ///
  /// @note This is intentionally not a regular setter. The particle id
  ///       is used to identify the whole particle. Setting it on an existing
  ///       particle is usually a mistake.
  Particle withParticleId(Barcode particleId) const {
    Particle p = *this;
    p.m_particleId = particleId;
    return p;
  }

  /// Set the process type that generated this particle.
  Particle &setProcess(ProcessType proc) { return m_process = proc, *this; }
  /// Set the space-time position four-vector.
  ///
  /// The component order is [x,y,z,t].
  Particle &setPosition4(const Vector4 &pos4) {
    m_position4 = pos4;
    return *this;
  }
  /// Set the space-time position four-vector from three-position and time.
  Particle &setPosition4(const Vector3 &position, Scalar time) {
    m_position4.head<3>() = position;
    m_position4[3] = time;
    return *this;
  }
  /// Set the space-time position four-vector from scalar components.
  Particle &setPosition4(Scalar x, Scalar y, Scalar z, Scalar time) {
    m_position4[0] = x;
    m_position4[1] = y;
    m_position4[2] = z;
    m_position4[3] = time;
    return *this;
  }
  /// Set the direction three-vector
  Particle &setDirection(const Vector3 &direction) {
    m_unitDirection = direction;
    m_unitDirection.normalize();
    return *this;
  }
  /// Set the direction three-vector from scalar components.
  Particle &setDirection(Scalar dx, Scalar dy, Scalar dz) {
    m_unitDirection[0] = dx;
    m_unitDirection[1] = dy;
    m_unitDirection[2] = dz;
    m_unitDirection.normalize();
    return *this;
  }
  /// Set the absolute momentum.
  Particle &setAbsMomentum(Scalar absMomentum) {
    m_absMomentum = absMomentum;
    return *this;
  }
  /// Change the energy by the given amount.
  ///
  /// Energy loss corresponds to a negative change. If the updated energy
  /// would result in an unphysical value, the particle is put to rest, i.e.
  /// its absolute momentum is set to zero.
  Particle &correctEnergy(Scalar delta) {
    const auto newEnergy = std::hypot(m_mass, m_absMomentum) + delta;
    if (newEnergy <= m_mass) {
      m_absMomentum = Scalar(0);
    } else {
      m_absMomentum = std::sqrt(newEnergy * newEnergy - m_mass * m_mass);
    }
    return *this;
  }

  /// Particle identifier within an event.
  constexpr Barcode particleId() const { return m_particleId; }
  /// Which type of process generated this particle.
  constexpr ProcessType process() const { return m_process; }
  /// PDG particle number that identifies the type.
  constexpr Acts::PdgParticle pdg() const { return m_pdg; }
  /// Particle charge.
  constexpr Scalar charge() const { return m_charge; }
  /// Particle mass.
  constexpr Scalar mass() const { return m_mass; }

  /// Space-time position four-vector.
  ///
  /// The component order is [x,y,z,t].
  constexpr const Vector4 &position4() const { return m_position4; }
  /// Three-position, i.e. spatial coordinates without the time.
  auto position() const { return m_position4.head<3>(); }
  /// Time coordinate.
  Scalar time() const { return m_position4[3]; }
  /// Energy-momentum four-vector.
  ///
  /// The component order is [px,py,pz,E].
  Vector4 momentum4() const {
    Vector4 mom4;
    // stored direction is always normalized
    mom4[0] = m_absMomentum * m_unitDirection[0];
    mom4[1] = m_absMomentum * m_unitDirection[1];
    mom4[2] = m_absMomentum * m_unitDirection[2];
    mom4[3] = energy();
    return mom4;
  }
  /// Unit three-direction, i.e. the normalized momentum three-vector.
  const Vector3 &unitDirection() const { return m_unitDirection; }
  /// Absolute momentum in the x-y plane.
  Scalar transverseMomentum() const {
    return m_absMomentum * m_unitDirection.head<2>().norm();
  }
  /// Absolute momentum.
  constexpr Scalar absMomentum() const { return m_absMomentum; }
  /// Total energy, i.e. norm of the four-momentum.
  Scalar energy() const { return std::hypot(m_mass, m_absMomentum); }

  /// Check if the particle is alive, i.e. is not at rest.
  constexpr operator bool() const { return Scalar(0) < m_absMomentum; }
  /// Check if the particle is dead, i.e is at rest.
  constexpr bool operator!() const { return m_absMomentum <= Scalar(0); }

  /// Set the material that the particle has passed.
  ///
  /// @param pathX0 passed material measured in radiation lengths
  /// @param pathL0 passed thickness measured in interaction lengths
  constexpr Particle &setMaterialPassed(Scalar pathX0, Scalar pathL0) {
    m_pathX0 = pathX0;
    m_pathL0 = pathL0;
    return *this;
  }
  /// Set the material limits.
  ///
  /// @param limitX0 maximum radiation lengths the particle can pass
  /// @param limitL0 maximum interaction lengths the particle can pass
  constexpr Particle &setMaterialLimits(Scalar limitX0, Scalar limitL0) {
    m_limitX0 = limitX0;
    m_limitL0 = limitL0;
    return *this;
  }
  /// The passed material measured in radiation lengths.
  constexpr Scalar pathInX0() const { return m_pathX0; }
  /// The passed material measured in interaction lengths.
  constexpr Scalar pathInL0() const { return m_pathL0; }
  /// The maximum radation length the particle is allowed to pass.
  constexpr Scalar pathLimitX0() const { return m_limitX0; }
  /// The maximum interaction length the particle is allowed to pass.
  constexpr Scalar pathLimitL0() const { return m_limitL0; }

 private:
  // identity, i.e. things that do not change over the particle lifetime.
  /// Particle identifier within the event.
  Barcode m_particleId;
  /// Process type specifier.
  ProcessType m_process = ProcessType::eUndefined;
  /// PDG particle number.
  Acts::PdgParticle m_pdg = Acts::PdgParticle::eInvalid;
  // Particle charge and mass.
  Scalar m_charge = Scalar(0);
  Scalar m_mass = Scalar(0);
  // kinematics, i.e. things that change over the particle lifetime.
  Vector3 m_unitDirection = Vector3::UnitZ();
  Scalar m_absMomentum = Scalar(0);
  Vector4 m_position4 = Vector4::Zero();
  // simulation-specific X0/L0 information and limits
  // these values are here to simplify the simulation of (nuclear) interactions.
  // instead of checking at every surface whether an interaction should occur we
  // can draw an overall limit once. the relevant interaction only needs to
  // be executed once the limit is reached.
  // this information is not really particle-specific and should probably be
  // handled separately. for now, storing it directly here is the simplest
  // solution.
  Scalar m_pathX0 = Scalar(0);
  Scalar m_pathL0 = Scalar(0);
  Scalar m_limitX0 = std::numeric_limits<Scalar>::max();
  Scalar m_limitL0 = std::numeric_limits<Scalar>::max();
};

std::ostream &operator<<(std::ostream &os, const Particle &particle);

}  // namespace ActsFatras
