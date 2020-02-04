// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

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
  /// @param id     Encoded identifier within an event
  /// @param pdg    PDG particle number
  /// @param mass   Particle mass in native units
  /// @param charge Particle charge in native units
  ///
  /// @warning It is the users responsibility that charge and mass match
  ///          the PDG particle number.
  Particle(Barcode id, Acts::PdgParticle pdg, Scalar mass, Scalar charge)
      : m_id(id), m_pdg(pdg), m_charge(charge), m_mass(mass) {}
  Particle(const Particle &) = default;
  Particle(Particle &&) = default;
  Particle &operator=(const Particle &) = default;
  Particle &operator=(Particle &&) = default;

  /// Set the space-time position four-vector.
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
    m_direction = direction;
    m_direction.normalize();
    return *this;
  }
  /// Set the direction three-vector from scalar components.
  Particle &setDirection(Scalar dx, Scalar dy, Scalar dz) {
    m_direction[0] = dx;
    m_direction[1] = dy;
    m_direction[2] = dz;
    m_direction.normalize();
    return *this;
  }
  /// Set the absolute momentum.
  Particle &setMomentum(Scalar momentum) {
    m_momentum = momentum;
    return *this;
  }
  /// Change the energy by the given amount.
  ///
  /// Energy loss corresponds to a negative change. If the updated energy
  /// would result in an unphysical value, the particle is put to rest, i.e.
  /// its absolute momentum is set to zero.
  Particle &correctEnergy(Scalar delta) {
    const auto newEnergy = std::hypot(m_mass, m_momentum) + delta;
    if (newEnergy <= m_mass) {
      m_momentum = Scalar(0);
    } else {
      m_momentum = std::sqrt(newEnergy * newEnergy - m_mass * m_mass);
    }
    return *this;
  }

  /// Encoded particle identifier within an event.
  Barcode id() const { return m_id; }
  /// PDG particle number that identifies the type.
  Acts::PdgParticle pdg() const { return m_pdg; }
  /// Particle charge.
  Scalar charge() const { return m_charge; }
  /// Particle mass.
  Scalar mass() const { return m_mass; }

  /// Space-time position four-vector.
  const Vector4 &position4() const { return m_position4; }
  /// Three-position, i.e. spatial coordinates without the time.
  auto position() const { return m_position4.head<3>(); }
  /// Time coordinate.
  Scalar time() const { return m_position4[3]; }
  /// Energy-momentum four-vector.
  Vector4 momentum4() const {
    Vector4 mom4;
    // stored direction is always normalized
    mom4[0] = m_momentum * m_direction[0];
    mom4[1] = m_momentum * m_direction[1];
    mom4[2] = m_momentum * m_direction[2];
    mom4[3] = energy();
    return mom4;
  }
  /// Three-direction, i.e. the normalized momentum three-vector.
  const Vector3 &direction() const { return m_direction; }
  /// Absolute momentum.
  Scalar momentum() const { return m_momentum; }
  /// Total energy, i.e. norm of the four-momentum.
  Scalar energy() const { return std::hypot(m_mass, m_momentum); }

  /// Charge over absolute momentum.
  Scalar chargeOverMomentum() const { return m_charge / m_momentum; }
  /// Relativistic velocity.
  Scalar beta() const { return m_momentum / energy(); }
  /// Relativistic gamma factor.
  Scalar gamma() const { return std::hypot(1, m_momentum / m_mass); }

  /// Check if the particle is alive, i.e. is not at rest.
  operator bool() const { return Scalar(0) < m_momentum; }
  /// Check if the particle is dead, i.e is at rest.
  bool operator!() const { return m_momentum <= Scalar(0); }

  /// Register material that the particle has passed.
  ///
  /// @param thicknessX0 material thickness measured in radiation lengths
  /// @param thicknessL0 material thickness measured in interaction lengths
  Particle &addPassedMaterial(Scalar thicknessX0, Scalar thicknessL0) {
    m_pathX0 += thicknessX0;
    m_pathL0 += thicknessL0;
    return *this;
  }
  /// Set the material limits.
  ///
  /// @param limitX0 maximum radiation lengths the particle can pass
  /// @param limitL0 maximum interaction lengths the particle can pass
  Particle &setMaterialLimits(Scalar limitX0, Scalar limitL0) {
    m_limitX0 = limitX0;
    m_limitL0 = limitL0;
    return *this;
  }
  /// The passed material measured in radiation lengths.
  Scalar pathInX0() const { return m_pathX0; }
  /// The passed material measured in interaction lengths.
  Scalar pathInL0() const { return m_pathL0; }
  /// The maximum radation length the particle is allowed to pass.
  Scalar pathLimitX0() const { return m_limitX0; }
  /// The maximum interaction length the particle is allowed to pass.
  Scalar pathLimitL0() const { return m_limitL0; }

 private:
  // identity, i.e. things that do not change over the particle lifetime.
  /// Particle identifier within the event.
  Barcode m_id;
  /// PDG particle number.
  Acts::PdgParticle m_pdg = Acts::PdgParticle::eInvalid;
  // Particle charge and mass.
  Scalar m_charge = Scalar(0);
  Scalar m_mass = Scalar(0);
  // kinematics, i.e. things that change over the particle lifetime.
  Vector3 m_direction = Vector3::UnitZ();
  Scalar m_momentum = Scalar(0);
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

}  // namespace ActsFatras
