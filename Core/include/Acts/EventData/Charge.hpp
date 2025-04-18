// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ChargeConcept.hpp"
#include "Acts/Utilities/Concepts.hpp"

#include <cassert>
#include <cmath>

namespace Acts {

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
  constexpr Neutral(float absQ) noexcept {
    assert((absQ == 0) && "Input charge must be zero");
    (void)absQ;
  }

  constexpr float absQ() const noexcept { return 0; }

  constexpr float extractCharge(ActsScalar /*qOverP*/) const noexcept {
    return 0.0f;
  }

  constexpr ActsScalar extractMomentum(ActsScalar qOverP) const noexcept {
    assert(qOverP >= 0 && "qOverP cannot be negative");
    return 1.0f / qOverP;
  }

  constexpr ActsScalar qOverP(ActsScalar momentum,
                              float signedQ) const noexcept {
    assert((signedQ != 0) && "charge must be 0");
    (void)signedQ;
    return 1.0f / momentum;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `Neutral` has no internal state.
  /// Must be available to provide a consistent interface.
  friend constexpr bool operator==(Neutral /*lhs*/, Neutral /*rhs*/) noexcept {
    return true;
  }
};

ACTS_STATIC_CHECK_CONCEPT(ChargeConcept, Neutral);

/// Charge and momentum interpretation for particles with +-e charge.
struct SinglyCharged {
  constexpr SinglyCharged() = default;

  // TODO remove this method after grad refactor; currently track parameters
  // depend on it
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  constexpr SinglyCharged(float absQ) noexcept {
    assert((absQ == UnitConstants::e) && "Input charge magnitude must be e");
    (void)absQ;
  }

  constexpr float absQ() const noexcept { return UnitConstants::e; }

  constexpr float extractCharge(ActsScalar qOverP) const noexcept {
    return std::copysign(UnitConstants::e, qOverP);
  }

  constexpr ActsScalar extractMomentum(ActsScalar qOverP) const noexcept {
    return extractCharge(qOverP) / qOverP;
  }

  constexpr ActsScalar qOverP(ActsScalar momentum,
                              float signedQ) const noexcept {
    assert((std::abs(signedQ) == UnitConstants::e) &&
           "absolute charge must be e");
    return signedQ / momentum;
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `SinglyCharged` has no internal state.
  /// Must be available to provide a consistent interface.
  friend constexpr bool operator==(SinglyCharged /*lhs*/,
                                   SinglyCharged /*rhs*/) noexcept {
    return true;
  }
};

ACTS_STATIC_CHECK_CONCEPT(ChargeConcept, SinglyCharged);

/// Charge and momentum interpretation for arbitrarily charged but not neutral
/// particles.
class NonNeutralCharge {
 public:
  /// Construct with the magnitude of the input charge.
  constexpr NonNeutralCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 < absQ) && "Input charge magnitude must be positive");
  }
  constexpr NonNeutralCharge(SinglyCharged /*unused*/) noexcept
      : m_absQ{UnitConstants::e} {}

  constexpr float absQ() const noexcept { return m_absQ; }

  constexpr float extractCharge(ActsScalar qOverP) const noexcept {
    return std::copysign(m_absQ, qOverP);
  }
  constexpr ActsScalar extractMomentum(ActsScalar qOverP) const noexcept {
    return extractCharge(qOverP) / qOverP;
  }

  constexpr ActsScalar qOverP(ActsScalar momentum,
                              float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absQ && "inconsistent charge");
    return signedQ / momentum;
  }

  /// Compare for equality.
  friend constexpr bool operator==(NonNeutralCharge lhs,
                                   NonNeutralCharge rhs) noexcept {
    return lhs.m_absQ == rhs.m_absQ;
  }

 private:
  float m_absQ{};
};

ACTS_STATIC_CHECK_CONCEPT(ChargeConcept, NonNeutralCharge);

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class AnyCharge {
 public:
  /// Construct with the magnitude of the input charge.
  constexpr AnyCharge(float absQ) noexcept : m_absQ{absQ} {
    assert((0 <= absQ) && "Input charge magnitude must be zero or positive");
  }
  constexpr AnyCharge(SinglyCharged /*unused*/) noexcept
      : m_absQ{UnitConstants::e} {}
  constexpr AnyCharge(Neutral /*unused*/) noexcept {}

  constexpr float absQ() const noexcept { return m_absQ; }

  constexpr float extractCharge(ActsScalar qOverP) const noexcept {
    return std::copysign(m_absQ, qOverP);
  }
  constexpr ActsScalar extractMomentum(ActsScalar qOverP) const noexcept {
    return (m_absQ != 0.0f) ? extractCharge(qOverP) / qOverP : 1.0f / qOverP;
  }

  constexpr ActsScalar qOverP(ActsScalar momentum,
                              float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absQ && "inconsistent charge");
    return (m_absQ != 0.0f) ? signedQ / momentum : 1.0f / momentum;
  }

  /// Compare for equality.
  friend constexpr bool operator==(AnyCharge lhs, AnyCharge rhs) noexcept {
    return lhs.m_absQ == rhs.m_absQ;
  }

 private:
  float m_absQ{};
};

ACTS_STATIC_CHECK_CONCEPT(ChargeConcept, AnyCharge);

/// @}

}  // namespace Acts
