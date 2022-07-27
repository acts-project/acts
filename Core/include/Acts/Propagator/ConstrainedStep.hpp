// This file is part of the Acts project.
//
// Copyright (C) 2018-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Acts {

/// A constrained step class for the steppers
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

  /// constructor from Scalar
  /// @param value is the user given initial value
  explicit ConstrainedStep(Scalar value)
      : m_direction(direction(value)) {
    m_values[user] = std::abs(value);
  }

  /// Assignment from one Scalar
  /// @note this will set only the accuracy, as this is the most
  /// exposed to the Propagator, this adapts also the direction
  ///
  /// @param value is the new accuracy value
  void setValue(Scalar value) {
    assert(value != not_set);
    /// set the accuracy value
    m_values[accuracy] = std::abs(value);
    // set/update the direction
    m_direction = direction(value);
    return (*this);
  }

  /// returns the min/max value depending on the direction
  Scalar value() const {
    Scalar value = (*std::min_element(m_values.begin(), m_values.end()));
    assert(value != not_set);
    return value * m_direction;
  }

  /// Access to a specific value
  ///
  /// @param type is the resquested parameter type
  Scalar value(Type type) const { return m_values[type] * m_direction; }

  Direction direction() const { return m_direction; }

  /// Access to currently leading min type
  ///
  Type currentType() const {
    return Type(std::min_element(m_values.begin(), m_values.end()) -
                m_values.begin());
  }

  /// release a certain constraint value
  /// to the (signed) biggest value available, hence
  /// it depends on the direction
  ///
  /// @param type is the constraint type to be released
  void release(Type type) {
    m_values[type] = not_set;
  }

  /// Update the step size of a certain type
  ///
  /// Only navigation and target abortion step size
  /// updates may change the sign due to overstepping
  ///
  /// @param value is the new value to be updated
  /// @param type is the constraint type
  /// @param releaseStep Allow step size to increase again
  void update(Scalar value, Type type, bool releaseStep = false) {
    assert(value != not_set);
    assert(direction(value) == m_direction);
    if (releaseStep) {
      release(type);
    }
    value = std::abs(value);
    // The check the current value and set it if appropriate
    Scalar cValue = m_values[type];
    m_values[type] = cValue < value ? cValue : value;
  }

  void scale(Scalar factor) {
    setValue(value() * factor);
  }

  /// return the split value as string for debugging
  std::string toString() const {
    std::stringstream dstream;

    // Helper method to avoid unreadable screen output
    auto streamValue = [&](Type type) {
      Scalar val = value(type);
      dstream << std::setw(5);
      if (std::abs(val) == not_set) {
        dstream << (val > 0 ? "+∞" : "-∞");
      } else {
        dstream << val;
      }
    };

    dstream << "(";
    streamValue(accuracy);
    dstream << ", ";
    streamValue(actor);
    dstream << ", ";
    streamValue(aborter);
    dstream << ", ";
    streamValue(user);
    dstream << ")";
    return dstream.str();
  }

private:
  inline static constexpr auto not_set = std::numeric_limits<Scalar>::max();

  /// the step size tuple
  std::array<Scalar, 4> m_values = {not_set, not_set, not_set, not_set};

  /// The Navigation direction
  NavigationDirection m_direction = NavigationDirection::Forward;
};

}  // namespace Acts
