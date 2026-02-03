// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ChargeConcept.hpp"

#include <cassert>
#include <cmath>

namespace Acts {

/// @ingroup eventdata
/// @defgroup eventdata-charge Charge interpretation for track parameters
///
/// Track parameters store a single coefficient that describes charge and
/// momentum. This is either charge/momentum or 1/momentum, but the
/// interpretation depends on what type of particle is described. In this code
/// base this coefficient is always referred to as `qOverP` (or
/// charge-over-momentum) even for uncharged particles. The following types are
/// used to restrict the particle charge magnitude (at compile time) and support
/// the umambigous extraction of charge and absolute momentum from said track
/// parameter coefficient.
///
/// All types are designed to be interchangeable. Each one can be
/// constructed with the input charge magnitude
///
/// ```cpp
/// Charge c(1_e);
/// ```
///
/// and can then be used to extract the charge value
///
/// ```cpp
/// auto q = c.extractCharge(qOverP);
/// ```
///
/// or the absolute momentum
///
/// ```cpp
/// auto p = c.extractMomentum(qOverP);
/// ```
///
/// from the charge-over-momentum track parameter.
///
/// @{

/// Charge and momentum interpretation for neutral particles.
struct Neutral {
  constexpr Neutral() = default;

  // TODO remove this method after grad refactor; currently track parameters
  // depend on it
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  /// @param absQ Absolute charge magnitude (must be zero for neutral particles)
  constexpr explicit Neutral(float absQ) noexcept {
    assert((absQ == 0) && "Input charge must be zero");
    static_cast<void>(absQ);
  }

  /// Get the absolute charge magnitude
  /// @return Always returns 0 for neutral particles
  constexpr float absQ() const noexcept { return 0; }

  /// Extract the signed charge from q/p
  /// @return Always returns 0 for neutral particles
  constexpr float extractCharge(double /*qOverP*/) const noexcept {
    return 0.0f;
  }

  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum (must be positive for neutral particles)
  /// @return Momentum magnitude calculated as 1/qOverP
  constexpr double extractMomentum(double qOverP) const noexcept {
    assert(qOverP >= 0 && "qOverP cannot be negative");
    return 1.0f / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must be 0 for neutral particles)
  /// @return Charge over momentum (1/momentum for neutral particles)
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert((signedQ != 0) && "charge must be 0");
    static_cast<void>(signedQ);
    return 1.0f / momentum;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `Neutral` has no internal state.
  /// Must be available to provide a consistent interface.
  constexpr bool operator==(const Neutral& rhs) const noexcept = default;
};

static_assert(ChargeConcept<Neutral>, "Neutral does not fulfill ChargeConcept");

/// Charge and momentum interpretation for particles with +-e charge.
struct SinglyCharged {
  constexpr SinglyCharged() = default;

  // TODO remove this method after grad refactor; currently track parameters
  // depend on it
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  /// @param absQ Absolute charge magnitude (must be e for singly charged particles)
  constexpr explicit SinglyCharged(float absQ) noexcept {
    assert((absQ == UnitConstants::e) && "Input charge magnitude must be e");
    static_cast<void>(absQ);
  }

  /// Get the absolute charge magnitude
  /// @return Elementary charge magnitude e
  constexpr float absQ() const noexcept { return UnitConstants::e; }

  /// Extract the signed charge from q/p
  /// @param qOverP Charge over momentum
  /// @return Signed elementary charge (+e or -e)
  constexpr float extractCharge(double qOverP) const noexcept {
    return static_cast<float>(std::copysign(UnitConstants::e, qOverP));
  }

  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum
  /// @return Momentum magnitude calculated as charge/qOverP
  constexpr double extractMomentum(double qOverP) const noexcept {
    return extractCharge(qOverP) / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must be ±e)
  /// @return Charge over momentum
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert((std::abs(signedQ) == UnitConstants::e) &&
           "absolute charge must be e");
    return signedQ / momentum;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `SinglyCharged` has no internal state.
  /// Must be available to provide a consistent interface.
  constexpr bool operator==(const SinglyCharged& rhs) const noexcept = default;
};

static_assert(ChargeConcept<SinglyCharged>,
              "SinglyCharged does not fulfill ChargeConcept");

/// Charge and momentum interpretation for arbitrarily charged but not neutral
/// particles.
class NonNeutralCharge {
 public:
  /// Construct with the magnitude of the input charge.
  /// @param absQ Absolute charge magnitude (must be positive for non-neutral particles)
  constexpr explicit NonNeutralCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 < absQ) && "Input charge magnitude must be positive");
  }
  /// Construct from a SinglyCharged particle
  constexpr explicit NonNeutralCharge(SinglyCharged /*unused*/) noexcept
      : m_absQ{UnitConstants::e} {}

  /// Get the absolute charge magnitude
  /// @return Absolute charge magnitude
  constexpr float absQ() const noexcept { return m_absQ; }

  /// Extract the signed charge from q/p
  /// @param qOverP Charge over momentum
  /// @return Signed charge with correct magnitude
  constexpr float extractCharge(double qOverP) const noexcept {
    return static_cast<float>(std::copysign(m_absQ, qOverP));
  }
  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum
  /// @return Momentum magnitude calculated as charge/qOverP
  constexpr double extractMomentum(double qOverP) const noexcept {
    return extractCharge(qOverP) / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must match stored charge magnitude)
  /// @return Charge over momentum
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absQ && "inconsistent charge");
    return signedQ / momentum;
  }

  /// Compare for equality.
  constexpr bool operator==(const NonNeutralCharge& rhs) const noexcept =
      default;

 private:
  float m_absQ{};
};

static_assert(ChargeConcept<NonNeutralCharge>,
              "NonNeutralCharge does not fulfill ChargeConcept");

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class AnyCharge {
 public:
  /// Construct with the magnitude of the input charge.
  /// @param absQ The absolute value of the charge magnitude
  constexpr explicit AnyCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 <= absQ) && "Input charge magnitude must be zero or positive");
  }
  /// Construct from a SinglyCharged particle
  constexpr explicit AnyCharge(SinglyCharged /*unused*/) noexcept
      : m_absQ{UnitConstants::e} {}
  /// Construct from a Neutral particle
  constexpr explicit AnyCharge(Neutral /*unused*/) noexcept {}

  /// Get the absolute charge magnitude
  /// @return Absolute charge magnitude (0 for neutral particles)
  constexpr float absQ() const noexcept { return m_absQ; }

  /// Extract the signed charge from q/p
  /// @param qOverP Charge over momentum
  /// @return Signed charge with correct magnitude (0 for neutral)
  constexpr float extractCharge(double qOverP) const noexcept {
    return static_cast<float>(std::copysign(m_absQ, qOverP));
  }
  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum
  /// @return Momentum magnitude (handles both charged and neutral particles)
  constexpr double extractMomentum(double qOverP) const noexcept {
    return (m_absQ != 0.0f) ? extractCharge(qOverP) / qOverP : 1.0 / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must match stored charge magnitude)
  /// @return Charge over momentum (handles both charged and neutral particles)
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absQ && "inconsistent charge");
    return (m_absQ != 0.0f) ? signedQ / momentum : 1.0 / momentum;
  }

  /// Compare for equality.
  constexpr bool operator==(const AnyCharge& rhs) const noexcept = default;

 private:
  float m_absQ{};
};

static_assert(ChargeConcept<AnyCharge>,
              "AnyCharge does not fulfill ChargeConcept");

/// @}

}  // namespace Acts
