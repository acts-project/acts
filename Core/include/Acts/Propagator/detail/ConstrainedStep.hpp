// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <array>
#include <iomanip>
#include <limits>
#include <sstream>

#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace detail {

/// A constrained step class for the steppers
struct ConstrainedStep {
  /// the types of constraints
  /// from accuracy - this can vary up and down given a good step estimator
  /// from actor    - this would be a typical navigation step
  /// from aborter  - this would be a target condition
  /// from user     - this is user given for what reason ever
  enum Type : int { accuracy = 0, actor = 1, aborter = 2, user = 3 };

  /// the step size tuple
  std::array<double, 4> values = {
      {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
       std::numeric_limits<double>::max(), std::numeric_limits<double>::max()}};

  /// The Navigation direction
  NavigationDirection direction = forward;

  /// Update the step size of a certain type
  ///
  /// Only navigation and target abortion step size
  /// updates may change the sign due to overstepping
  ///
  /// @param value is the new value to be updated
  /// @param type is the constraint type
  void update(const double& value, Type type, bool releaseStep = false) {
    if (releaseStep) {
      release(type);
    }
    // The check the current value and set it if appropriate
    double cValue = values[type];
    values[type] = cValue * cValue < value * value ? cValue : value;
  }

  /// release a certain constraint value
  /// to the (signed) biggest value available, hence
  /// it depends on the direction
  ///
  /// @param type is the constraint type to be released
  void release(Type type) {
    double mvalue = (direction == forward)
                        ? (*std::max_element(values.begin(), values.end()))
                        : (*std::min_element(values.begin(), values.end()));
    values[type] = mvalue;
  }

  /// constructor from double
  /// @paramn value is the initial value for all steps
  ConstrainedStep(double value)
      : values({{value, value, value, value}}),
        direction(value > 0. ? forward : backward) {}

  /// The assignment operator from one double
  /// @note this will set only the accuracy, as this is the most
  /// exposed to the Propagator, this adapts also the direction
  ///
  /// @param value is the new accuracy value
  ConstrainedStep& operator=(const double& value) {
    /// set the accuracy value
    values[accuracy] = value;
    // set/update the direction
    direction = value > 0. ? forward : backward;
    return (*this);
  }

  /// Cast operator to double, returning the min/max value
  /// depending on the direction
  operator double() const {
    if (direction == forward) {
      return (*std::min_element(values.begin(), values.end()));
    }
    return (*std::max_element(values.begin(), values.end()));
  }

  /// Access to a specific value
  ///
  /// @param type is the resquested parameter type
  double value(Type type) const { return values[type]; }

  /// Return the maximum step constraint
  /// @return The max step constraint
  double max() const {
    return (*std::max_element(values.begin(), values.end()));
  }

  /// Return the minimum step constraint
  /// @return The min step constraint
  double min() const {
    return (*std::min_element(values.begin(), values.end()));
  }

  /// Access to currently leading min type
  ///
  Type currentType() const {
    if (direction == forward) {
      return Type(std::min_element(values.begin(), values.end()) -
                  values.begin());
    }
    return Type(std::max_element(values.begin(), values.end()) -
                values.begin());
  }

  /// return the split value as string for debugging
  std::string toString() const;
};

inline std::string ConstrainedStep::toString() const {
  std::stringstream dstream;
  dstream << "(" << std::setw(5) << values[accuracy];
  dstream << ", " << std::setw(5) << values[actor];
  dstream << ", " << std::setw(5) << values[aborter];
  dstream << ", " << std::setw(5) << values[user] << " )";
  return dstream.str();
}

}  // namespace detail
}  // namespace Acts
