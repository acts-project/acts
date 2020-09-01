// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace Acts {
// This is the custom error code enum
enum class PropagatorError {
  Failure = 1,
  WrongDirection = 2,
  StepCountLimitReached = 3
};

namespace detail {
// Define a custom error code category derived from std::error_category
class PropagatorErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "PropagatorError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<PropagatorError>(c)) {
      case PropagatorError::Failure:
        return "Propagation failed";
      case PropagatorError::WrongDirection:
        return "Propagation occurred in the wrong direction";
      case PropagatorError::StepCountLimitReached:
        return "Propagation reached the configured maximum number of steps";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::PropagatorErrorCategory& PropagatorErrorCategory() {
  static detail::PropagatorErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::PropagatorError e) {
  return {static_cast<int>(e), Acts::PropagatorErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::PropagatorError> : std::true_type {};
}  // namespace std
