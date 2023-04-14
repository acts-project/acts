// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Direction.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>

namespace Acts {

/// A constrained step class for the steppers
///
/// As simple as this class looks it hides a few very important details:
/// - Overstepping handling. The step size sign will flip if we happened to pass
/// our target.
/// - Convergence handling. Smaller and smaller step sizes have to be used in
/// order to converge on a target.
///
/// Because of the points mentioned above, the update function will always
/// prefer step sizes that point opposite the nagivation direction. A side
/// effect of this is that we will propagate in the opposite direction if the
/// target is "behind us".
///
/// The hierarchy is:
/// - Overstepping resolution / backpropagation
/// - Convergence
/// - Step into the void with `std::numeric_limits<Scalar>::max()`
class ConstrainedStep {
 public:
  using Scalar = ActsScalar;

  /// the types of constraints
  /// from accuracy - this can vary up and down given a good step estimator
  /// from actor    - this would be a typical navigation step
  /// from aborter  - this would be a target condition
  /// from user     - this is user given for what reason ever
  enum Type : int { accuracy = 0, actor = 1, aborter = 2, user = 3 };

  /// Number of iterations needed by the stepsize finder
  /// (e.g. Runge-Kutta) of the stepper.
  size_t nStepTrials = std::numeric_limits<size_t>::max();

  constexpr ConstrainedStep() = default;

  /// constructor from Scalar
  /// navigation direction is inferred by the sign of the step size
  /// @param value is the user given initial value
  constexpr explicit ConstrainedStep(Scalar value) {
    m_values[user] = std::abs(value);
    m_direction = Direction::fromScalar(value);
  }

  /// set accuracy by one Scalar
  ///
  /// this will set only the accuracy, as this is the most
  /// exposed to the Propagator
  ///
  /// @param value is the new accuracy value
  constexpr void setValue(Scalar value) {
    /// set the accuracy value
    m_values[accuracy] = value * m_direction;
  }

  /// returns the min step size
  constexpr Scalar value() const { return value(currentType()); }

  /// Access a specific value
  ///
  /// @param type is the requested parameter type
  constexpr Scalar value(Type type) const {
    return m_values[type] * m_direction;
  }

  /// Access the currently leading type
  constexpr Type currentType() const {
    return Type(std::min_element(m_values.begin(), m_values.end()) -
                m_values.begin());
  }

  /// release a certain constraint value
  ///
  /// @param type is the constraint type to be released
  constexpr void release(Type type) { m_values[type] = kNotSet; }

  /// Update the step size of a certain type
  ///
  /// Only navigation and target abortion step size
  /// updates may change the sign due to overstepping
  ///
  /// @param value is the new value to be updated
  /// @param type is the constraint type
  /// @param releaseStep Allow step size to increase again
  constexpr void update(Scalar value, Type type, bool releaseStep = false) {
    if (releaseStep) {
      release(type);
    }
    // check the current value and set it if appropriate
    // this will also allow signed values due to overstepping
    if (std::abs(value) <= std::abs(m_values[type])) {
      m_values[type] = value * m_direction;
    }
  }

  constexpr void scale(Scalar factor) {
    assert(factor > 0 && "ConstrainedStep scale factor was zero or negative.");
    m_values[accuracy] = value() * factor * m_direction;
  }

  std::ostream& toStream(std::ostream& os) const {
    // Helper method to avoid unreadable screen output
    auto streamValue = [&](Type type) {
      Scalar val = value(type);
      os << std::setw(5);
      if (std::abs(val) == kNotSet) {
        os << (val > 0 ? "+∞" : "-∞");
      } else {
        os << val;
      }
    };

    os << "(";
    streamValue(accuracy);
    os << ", ";
    streamValue(actor);
    os << ", ";
    streamValue(aborter);
    os << ", ";
    streamValue(user);
    os << ")";

    return os;
  }

  std::string toString() const {
    std::stringstream dstream;
    toStream(dstream);
    return dstream.str();
  }

 private:
  inline static constexpr auto kNotSet = std::numeric_limits<Scalar>::max();

  /// the step size tuple
  /// all values point in the `m_direction`
  std::array<Scalar, 4> m_values = {kNotSet, kNotSet, kNotSet, kNotSet};
  /// the navigation direction
  /// the direction is invariant after initialization
  Direction m_direction = Direction::Forward;
};

inline std::ostream& operator<<(std::ostream& os, const ConstrainedStep& step) {
  return step.toStream(os);
}

}  // namespace Acts
