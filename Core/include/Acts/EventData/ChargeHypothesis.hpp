// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cassert>
#include <cmath>

namespace Acts {

/// @ingroup eventdata
/// @defgroup eventdata-charge Charge hypothesis for track reconstruction
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
/// ChargeHypothesis c(1_e);
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

/// Charge and momentum interpretation for arbitrarily charged particles.
///
/// Only a charge magnitude identical to zero is interpreted as representing a
/// neutral particle. This avoids ambiguities that might arise from using an
/// approximate comparison with an arbitrary epsilon.
class ChargeHypothesis final {
 public:
  /// Construct with the magnitude of the input charge.
  /// @param absoluteCharge The absolute value of the charge magnitude
  constexpr explicit ChargeHypothesis(float absoluteCharge) noexcept
      : m_absoluteCharge{absoluteCharge} {
    assert(absoluteCharge >= 0 &&
           "Input charge magnitude must be zero or positive");
  }

  /// Get the absolute charge magnitude
  /// @return Absolute charge magnitude (0 for neutral particles)
  constexpr float absoluteCharge() const noexcept { return m_absoluteCharge; }

  /// Extract the signed charge from q/p
  /// @param qOverP Charge over momentum
  /// @return Signed charge with correct magnitude (0 for neutral)
  constexpr float extractCharge(double qOverP) const noexcept {
    return static_cast<float>(std::copysign(m_absoluteCharge, qOverP));
  }
  /// Extract momentum magnitude from q/p
  /// @param qOverP Charge over momentum
  /// @return Momentum magnitude (handles both charged and neutral particles)
  constexpr double extractMomentum(double qOverP) const noexcept {
    return (m_absoluteCharge != 0.0f) ? extractCharge(qOverP) / qOverP
                                      : 1.0 / qOverP;
  }

  /// Compute q/p from momentum and signed charge
  /// @param momentum Particle momentum magnitude
  /// @param signedQ Signed charge (must match stored charge magnitude)
  /// @return Charge over momentum (handles both charged and neutral particles)
  constexpr double qOverP(double momentum, float signedQ) const noexcept {
    assert(std::abs(signedQ) == m_absoluteCharge && "inconsistent charge");
    return (m_absoluteCharge != 0.0f) ? signedQ / momentum : 1.0 / momentum;
  }

  /// Compare for equality.
  /// @param rhs The ChargeHypothesis to compare to
  /// @return True if the two ChargeHypothesis objects are equal, false otherwise
  constexpr bool operator==(const ChargeHypothesis& rhs) const noexcept =
      default;

 private:
  float m_absoluteCharge{};
};

/// @}

}  // namespace Acts
