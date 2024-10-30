// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ParticleOutcome.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <cmath>
#include <iosfwd>
#include <optional>

namespace ActsFatras {

/// Particle identity information and kinematic state.
///
/// Also stores some simulation-specific properties.
class Particle {
 public:
  using Scalar = Acts::ActsScalar;
  using Vector3 = Acts::ActsVector<3>;
  using Vector4 = Acts::ActsVector<4>;

  struct State {
    // kinematics, i.e. things that change over the particle lifetime.
    Vector4 position4 = Vector4::Zero();
    Vector3 direction = Vector3::Zero();
    Scalar absMomentum = Scalar{0};
    /// proper time in the particle rest frame
    Scalar properTime = Scalar{0};
    // accumulated material
    Scalar pathInX0 = Scalar{0};
    Scalar pathInL0 = Scalar{0};
    /// number of hits
    std::uint32_t numberOfHits = 0;
    /// reference surface
    const Acts::Surface *referenceSurface{nullptr};
  };

  template <bool ReadOnly>
  class StateProxyImpl {
   public:
    using ParticleType = std::conditional_t<ReadOnly, const Particle, Particle>;
    using StateType = std::conditional_t<ReadOnly, const State, State>;

    StateProxyImpl(ParticleType &particle, StateType &state)
        : m_particle(&particle), m_state(&state) {}

    const Particle &particle() const { return *m_particle; }
    Particle &particle()
      requires(!ReadOnly)
    {
      return *m_particle;
    }

    const Vector4 &fourPosition() const { return m_state->position4; }
    const Vector3 &direction() const { return m_state->direction; }
    Scalar absoluteMomentum() const { return m_state->absMomentum; }
    Scalar properTime() const { return m_state->properTime; }
    Scalar pathInX0() const { return m_state->pathInX0; }
    Scalar pathInL0() const { return m_state->pathInL0; }
    std::uint32_t numberOfHits() const { return m_state->numberOfHits; }
    const Acts::Surface *referenceSurface() const {
      return m_state->referenceSurface;
    }

    Vector4 &fourPosition()
      requires(!ReadOnly)
    {
      return m_state->position4;
    }
    Vector3 &direction()
      requires(!ReadOnly)
    {
      return m_state->direction;
    }
    Scalar &absMomentum()
      requires(!ReadOnly)
    {
      return m_state->absMomentum;
    }
    Scalar &properTime()
      requires(!ReadOnly)
    {
      return m_state->properTime;
    }
    Scalar &pathInX0()
      requires(!ReadOnly)
    {
      return m_state->pathInX0;
    }
    Scalar &pathInL0()
      requires(!ReadOnly)
    {
      return m_state->pathInL0;
    }
    std::uint32_t &numberOfHits()
      requires(!ReadOnly)
    {
      return m_state->numberOfHits;
    }
    const Acts::Surface *&referenceSurface()
      requires(!ReadOnly)
    {
      return m_state->referenceSurface;
    }

    template <bool OtherReadOnly>
    void copyFrom(const StateProxyImpl<OtherReadOnly> &other)
      requires(!ReadOnly)
    {
      m_state->position4 = other.fourPosition();
      m_state->direction = other.direction();
      m_state->absMomentum = other.absoluteMomentum();
      m_state->properTime = other.properTime();
      m_state->pathInX0 = other.pathInX0();
      m_state->pathInL0 = other.pathInL0();
      m_state->numberOfHits = other.numberOfHits();
      m_state->referenceSurface = other.referenceSurface();
    }

    /// Check if the particle has a reference surface.
    bool hasReferenceSurface() const {
      return m_state->referenceSurface != nullptr;
    }

    /// Three-position, i.e. spatial coordinates without the time.
    Vector3 position() const {
      return m_state->position4.template segment<3>(Acts::ePos0);
    }
    /// Time coordinate.
    Scalar time() const { return m_state->position4[Acts::eTime]; }
    /// Energy-momentum four-vector.
    Vector4 fourMomentum() const {
      Vector4 mom4;
      // stored direction is always normalized
      mom4[Acts::eMom0] =
          m_state->absMomentum * m_state->direction[Acts::ePos0];
      mom4[Acts::eMom1] =
          m_state->absMomentum * m_state->direction[Acts::ePos1];
      mom4[Acts::eMom2] =
          m_state->absMomentum * m_state->direction[Acts::ePos2];
      mom4[Acts::eEnergy] = energy();
      return mom4;
    }
    /// Polar angle.
    Scalar theta() const {
      return Acts::VectorHelpers::theta(m_state->direction);
    }
    /// Azimuthal angle.
    Scalar phi() const { return Acts::VectorHelpers::phi(m_state->direction); }
    /// Absolute momentum in the x-y plane.
    Scalar transverseMomentum() const {
      return m_state->absMomentum * Acts::VectorHelpers::perp(direction());
    }
    /// Momentum vector.
    Vector3 momentum() const {
      return m_state->absMomentum * m_state->direction;
    }
    /// Total energy, i.e. norm of the four-momentum.
    Scalar energy() const {
      return Acts::fastHypot(m_particle->mass(), m_state->absMomentum);
    }
    /// Particl qOverP.
    Scalar qOverP() const {
      return m_particle->hypothesis().qOverP(absoluteMomentum(),
                                             m_particle->charge());
    }
    /// Eta.
    Scalar eta() const { return Acts::VectorHelpers::eta(direction()); }

    /// Set the space-time position four-vector.
    StateProxyImpl<false> &setPosition4(const Vector4 &pos4)
      requires(!ReadOnly)
    {
      m_state->position4 = pos4;
      return *this;
    }
    /// Set the space-time position four-vector from three-position and time.
    StateProxyImpl<false> &setPosition4(const Vector3 &position, Scalar time)
      requires(!ReadOnly)
    {
      m_state->position4.template segment<3>(Acts::ePos0) = position;
      m_state->position4[Acts::eTime] = time;
      return *this;
    }
    /// Set the space-time position four-vector from scalar components.
    StateProxyImpl<false> &setPosition4(Scalar x, Scalar y, Scalar z,
                                        Scalar time)
      requires(!ReadOnly)
    {
      m_state->position4[Acts::ePos0] = x;
      m_state->position4[Acts::ePos1] = y;
      m_state->position4[Acts::ePos2] = z;
      m_state->position4[Acts::eTime] = time;
      return *this;
    }
    /// Set the direction three-vector
    StateProxyImpl<false> &setDirection(const Vector3 &direction)
      requires(!ReadOnly)
    {
      m_state->direction = direction;
      m_state->direction.normalize();
      return *this;
    }
    /// Set the direction three-vector from scalar components.
    StateProxyImpl<false> &setDirection(Scalar dx, Scalar dy, Scalar dz)
      requires(!ReadOnly)
    {
      m_state->direction[Acts::ePos0] = dx;
      m_state->direction[Acts::ePos1] = dy;
      m_state->direction[Acts::ePos2] = dz;
      m_state->direction.normalize();
      return *this;
    }
    /// Set the absolute momentum.
    StateProxyImpl<false> &setAbsoluteMomentum(Scalar absMomentum)
      requires(!ReadOnly)
    {
      m_state->absMomentum = absMomentum;
      return *this;
    }
    /// Set the accumulated material measured in radiation/interaction lengths.
    ///
    /// @param pathInX0 accumulated material measured in radiation lengths
    /// @param pathInL0 accumulated material measured in interaction lengths
    StateProxyImpl<false> &setMaterialPassed(Scalar pathInX0, Scalar pathInL0)
      requires(!ReadOnly)
    {
      m_state->pathInX0 = pathInX0;
      m_state->pathInL0 = pathInL0;
      return *this;
    }
    /// Set the proper time in the particle rest frame.
    ///
    /// @param properTime passed proper time in the rest frame
    StateProxyImpl<false> &setProperTime(Scalar properTime)
      requires(!ReadOnly)
    {
      m_state->properTime = properTime;
      return *this;
    }
    /// Change the energy by the given amount.
    ///
    /// Energy loss corresponds to a negative change. If the updated energy
    /// would result in an unphysical value, the particle is put to rest, i.e.
    /// its absolute momentum is set to zero.
    StateProxyImpl<false> &correctEnergy(Scalar delta)
      requires(!ReadOnly)
    {
      const auto newEnergy =
          Acts::fastHypot(m_particle->mass(), m_state->absMomentum) + delta;
      if (newEnergy <= m_particle->mass()) {
        m_state->absMomentum = Scalar{0};
      } else {
        m_state->absMomentum = std::sqrt(
            newEnergy * newEnergy - m_particle->mass() * m_particle->mass());
      }
      return *this;
    }
    StateProxyImpl<false> &setNumberOfHits(std::uint32_t numberOfHits)
      requires(!ReadOnly)
    {
      m_state->numberOfHits = numberOfHits;
      return *this;
    }
    /// Set the reference surface.
    ///
    /// @param surface reference surface
    StateProxyImpl<false> &setReferenceSurface(const Acts::Surface *surface)
      requires(!ReadOnly)
    {
      m_state->referenceSurface = surface;
      return *this;
    }

    /// Bound track parameters.
    Acts::Result<Acts::BoundTrackParameters> boundParameters(
        const Acts::GeometryContext &gctx) const {
      if (!m_particle->hasReferenceSurface()) {
        return Acts::Result<Acts::BoundTrackParameters>::failure(
            std::error_code());
      }
      Acts::Result<Acts::Vector2> localResult =
          m_particle->m_referenceSurface->globalToLocal(gctx, position(),
                                                        direction());
      if (!localResult.ok()) {
        return localResult.error();
      }
      Acts::BoundVector params;
      params << localResult.value(), phi(), theta(), qOverP(), time();
      return Acts::BoundTrackParameters(
          m_particle->referenceSurface()->getSharedPtr(), params, std::nullopt,
          m_particle->hypothesis());
    }

    /// Curvilinear track parameters.
    Acts::CurvilinearTrackParameters curvilinearParameters() const {
      return Acts::CurvilinearTrackParameters(fourPosition(), direction(),
                                              qOverP(), std::nullopt,
                                              m_particle->hypothesis());
    }

   private:
    ParticleType *m_particle;
    StateType *m_state;
  };
  using StateProxy = StateProxyImpl<false>;
  using ConstStateProxy = StateProxyImpl<true>;

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
  Particle &setProcess(ProcessType proc) {
    m_process = proc;
    return *this;
  }
  /// Set the pdg.
  Particle setPdg(Acts::PdgParticle pdg) {
    m_pdg = pdg;
    return *this;
  }
  /// Set the charge.
  Particle setCharge(Scalar charge) {
    m_charge = charge;
    return *this;
  }
  /// Set the mass.
  Particle setMass(Scalar mass) {
    m_mass = mass;
    return *this;
  }
  /// Set the particle ID.
  Particle &setParticleId(Barcode barcode) {
    m_particleId = barcode;
    return *this;
  }

  /// Particle identifier within an event.
  Barcode particleId() const { return m_particleId; }
  /// Which type of process generated this particle.
  ProcessType process() const { return m_process; }
  /// PDG particle number that identifies the type.
  Acts::PdgParticle pdg() const { return m_pdg; }
  /// Absolute PDG particle number that identifies the type.
  Acts::PdgParticle absolutePdg() const {
    return Acts::makeAbsolutePdgParticle(pdg());
  }
  /// Particle charge.
  Scalar charge() const { return m_charge; }
  /// Particle absolute charge.
  Scalar absoluteCharge() const { return std::abs(m_charge); }
  /// Particle mass.
  Scalar mass() const { return m_mass; }

  /// Particle hypothesis.
  Acts::ParticleHypothesis hypothesis() const {
    return Acts::ParticleHypothesis(absolutePdg(), mass(), absoluteCharge());
  }

  bool isSecondary() const {
    return particleId().vertexSecondary() != 0 ||
           particleId().generation() != 0 || particleId().subParticle() != 0;
  }

  /// Set the outcome of particle.
  ///
  /// @param outcome outcome code
  Particle &setOutcome(ParticleOutcome outcome) {
    m_outcome = outcome;
    return *this;
  }

  /// Particle outcome.
  ParticleOutcome outcome() const { return m_outcome; }

  ConstStateProxy initialState() const {
    return ConstStateProxy(*this, m_initialState);
  }
  ConstStateProxy lastState() const {
    return ConstStateProxy(*this, m_lastState);
  }

  StateProxy initialState() { return StateProxy(*this, m_initialState); }
  StateProxy lastState() { return StateProxy(*this, m_lastState); }

  template <typename F>
  Particle &initialState(F &&f) {
    f(initialState());
    return *this;
  }
  template <typename F>
  Particle &lastState(F &&f) {
    f(lastState());
    return *this;
  }

 private:
  // identity, i.e. things that do not change over the particle lifetime.
  /// Particle identifier within the event.
  Barcode m_particleId;
  /// Process type specifier.
  ProcessType m_process = ProcessType::eUndefined;
  /// PDG particle number.
  Acts::PdgParticle m_pdg = Acts::PdgParticle::eInvalid;
  // Particle charge and mass.
  Scalar m_charge = Scalar{0};
  Scalar m_mass = Scalar{0};

  State m_initialState;
  State m_lastState;

  // additional final information
  /// outcome
  ParticleOutcome m_outcome = ParticleOutcome::Unknown;
};

std::ostream &operator<<(std::ostream &os, const Particle &particle);

}  // namespace ActsFatras
