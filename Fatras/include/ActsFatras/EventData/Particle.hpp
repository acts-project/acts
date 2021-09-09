// This file is part of the Acts project.
//
// Copyright (C) 2018-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Utilities/PdgParticle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cmath>
#include <iosfwd>
#include <limits>

namespace ActsFatras {

/// Particle identity information and kinematic state.
///
/// Also stores some simulation-specific properties.
class Particle {
 public:
  using Scalar = Acts::ActsScalar;
  using Vector3 = Acts::ActsVector<3>;
  using Vector4 = Acts::ActsVector<4>;

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
  /// Set the particle ID.
  Particle &setParticleId(Barcode barcode) {
    return m_particleId = barcode, *this;
  }
  /// Set the space-time position four-vector.
  Particle &setPosition4(const Vector4 &pos4) {
    m_position4 = pos4;
    return *this;
  }
  /// Set the space-time position four-vector from three-position and time.
  Particle &setPosition4(const Vector3 &position, Scalar time) {
    m_position4.segment<3>(Acts::ePos0) = position;
    m_position4[Acts::eTime] = time;
    return *this;
  }
  /// Set the space-time position four-vector from scalar components.
  Particle &setPosition4(Scalar x, Scalar y, Scalar z, Scalar time) {
    m_position4[Acts::ePos0] = x;
    m_position4[Acts::ePos1] = y;
    m_position4[Acts::ePos2] = z;
    m_position4[Acts::eTime] = time;
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
    m_unitDirection[Acts::ePos0] = dx;
    m_unitDirection[Acts::ePos1] = dy;
    m_unitDirection[Acts::ePos2] = dz;
    m_unitDirection.normalize();
    return *this;
  }
  /// Set the absolute momentum.
  Particle &setAbsoluteMomentum(Scalar absMomentum) {
    m_absMomentum = absMomentum;
    return *this;
  }

  /// Set the particle charge.
  Particle &setCharge(Scalar charge) { return m_charge = charge, *this; }

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
  /// PDG particle number of its parent.
  constexpr Acts::PdgParticle parentPDG() const { return m_parent_pdg; }
  Particle &setParentPDG(Acts::PdgParticle ppid) { return m_parent_pdg = ppid, *this; }

  /// Space-time position four-vector.
  constexpr const Vector4 &fourPosition() const { return m_position4; }
  /// Three-position, i.e. spatial coordinates without the time.
  auto position() const { return m_position4.segment<3>(Acts::ePos0); }
  /// Time coordinate.
  Scalar time() const { return m_position4[Acts::eTime]; }
  /// Energy-momentum four-vector.
  Vector4 fourMomentum() const {
    Vector4 mom4;
    // stored direction is always normalized
    mom4[Acts::eMom0] = m_absMomentum * m_unitDirection[Acts::ePos0];
    mom4[Acts::eMom1] = m_absMomentum * m_unitDirection[Acts::ePos1];
    mom4[Acts::eMom2] = m_absMomentum * m_unitDirection[Acts::ePos2];
    mom4[Acts::eEnergy] = energy();
    return mom4;
  }
  /// Unit three-direction, i.e. the normalized momentum three-vector.
  const Vector3 &unitDirection() const { return m_unitDirection; }
  /// Absolute momentum in the x-y plane.
  Scalar transverseMomentum() const {
    return m_absMomentum * m_unitDirection.segment<2>(Acts::eMom0).norm();
  }
  /// Absolute momentum.
  constexpr Scalar absoluteMomentum() const { return m_absMomentum; }
  /// Total energy, i.e. norm of the four-momentum.
  Scalar energy() const { return std::hypot(m_mass, m_absMomentum); }

  /// Check if the particle is alive, i.e. is not at rest.
  constexpr operator bool() const { return Scalar(0) < m_absMomentum; }
  /// Check if the particle is dead, i.e is at rest.
  constexpr bool operator!() const { return m_absMomentum <= Scalar(0); }

  // simulation specific properties

  /// Set the proper time in the particle rest frame.
  ///
  /// @param properTime passed proper time in the rest frame
  constexpr Particle &setProperTime(Scalar properTime) {
    m_properTime = properTime;
    return *this;
  }
  /// Proper time in the particle rest frame.
  constexpr Scalar properTime() const { return m_properTime; }

  /// Set the accumulated material measured in radiation/interaction lengths.
  ///
  /// @param pathInX0 accumulated material measured in radiation lengths
  /// @param pathInL0 accumulated material measured in interaction lengths
  constexpr Particle &setMaterialPassed(Scalar pathInX0, Scalar pathInL0) {
    m_pathInX0 = pathInX0;
    m_pathInL0 = pathInL0;
    return *this;
  }
  /// Accumulated path within material measured in radiation lengths.
  constexpr Scalar pathInX0() const { return m_pathInX0; }
  /// Accumulated path within material measured in interaction lengths.
  constexpr Scalar pathInL0() const { return m_pathInL0; }

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
  // proper time in the particle rest frame
  Scalar m_properTime = Scalar(0);
  // accumulated material
  Scalar m_pathInX0 = Scalar(0);
  Scalar m_pathInL0 = Scalar(0);

  // PDG particle number for its parent
  Acts::PdgParticle m_parent_pdg = Acts::PdgParticle::eInvalid;
};

std::ostream &operator<<(std::ostream &os, const Particle &particle);

}  // namespace ActsFatras
