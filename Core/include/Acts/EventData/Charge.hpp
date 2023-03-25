// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"

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
  Neutral() = default;
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  template <typename T>
  constexpr Neutral(T absQ) noexcept {
    assert((absQ == static_cast<T>(0)) and "Input charge must be zero");
    // suppress `unused variable` warning in non-debug builds
    (void)(absQ);
  }

  template <typename T>
  constexpr T extractCharge(T /* pInv */) const noexcept {
    return 0;
  }
  template <typename T>
  constexpr T extractMomentum(T pInv) const noexcept {
    // the abs is not strictly needed. it is added to protect against invalid,
    // i.e. negative, 1/p values to ensure that the output is still correct.
    return std::abs(1 / pInv);
  }

  /// Compare for equality.
  ///
  /// This is always `true` as `Neutral` has no internal state.
  /// Must be available to provide a consistent interface.
  friend constexpr bool operator==(Neutral /*lhs*/, Neutral /*rhs*/) noexcept {
    return true;
  }
};

/// Charge and momentum interpretation for particles with +-e charge.
struct SinglyCharged {
  SinglyCharged() = default;
  /// Construct and verify the input charge magnitude (in debug builds).
  ///
  /// This constructor is only provided to allow consistent construction.
  template <typename T>
  constexpr SinglyCharged(T absQ) noexcept {
    assert((absQ == static_cast<T>(UnitConstants::e)) and
           "Input charge magnitude must be e");
    // suppress `unused variable` warning in non-debug builds
    (void)(absQ);
  }

  template <typename T>
  constexpr T extractCharge(T qOverP) const noexcept {
    return std::copysign(static_cast<T>(UnitConstants::e), qOverP);
  }
  template <typename T>
  constexpr T extractMomentum(T qOverP) const noexcept {
    return std::abs(static_cast<T>(UnitConstants::e) / qOverP);
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

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class AnyCharge {
 public:
  /// Delete default constructor to ensure charge is always explicitely given.
  AnyCharge() = delete;
  /// Construct with the magnitude of the input charge.
  template <typename T>
  constexpr AnyCharge(T absQ) noexcept : m_magnitude(std::abs(absQ)) {
    assert((0 <= absQ) and "Input charge magnitude must be zero or positive");
  }

  template <typename T>
  constexpr T extractCharge(T qOverP) const noexcept {
    return std::copysign(static_cast<T>(m_magnitude), qOverP);
  }
  template <typename T>
  constexpr T extractMomentum(T qOverP) const noexcept {
    return (m_magnitude != 0.0f)
               ? std::abs(static_cast<T>(m_magnitude) / qOverP)
               : std::abs(1 / qOverP);
  }

 private:
  float m_magnitude;

  /// Compare for equality.
  friend constexpr bool operator==(AnyCharge lhs, AnyCharge rhs) noexcept {
    return lhs.m_magnitude == rhs.m_magnitude;
  }
};

/// @}

}  // namespace Acts
