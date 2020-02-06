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
enum class EigenStepperError {
  StepSizeStalled = 1,
  StepInvalid = 2,
  StepSizeAdjustmentFailed = 3
};

namespace detail {
// Define a custom error code category derived from std::error_category
class EigenStepperErrorCategory : public std::error_category {
 public:
  // Return a short descriptive name for the category
  const char* name() const noexcept final { return "EigenStepperError"; }
  // Return what each enum means in text
  std::string message(int c) const final {
    switch (static_cast<EigenStepperError>(c)) {
      case EigenStepperError::StepSizeStalled:
        return "Step size fell below minimum threshold";
      case EigenStepperError::StepInvalid:
        return "Step calculation was invalid";
      case EigenStepperError::StepSizeAdjustmentFailed:
        return "Step size adjustment exceeds maximum trials";
      default:
        return "unknown";
    }
  }
};
}  // namespace detail

// Declare a global function returning a static instance of the custom category
extern inline const detail::EigenStepperErrorCategory&
EigenStepperErrorCategory() {
  static detail::EigenStepperErrorCategory c;
  return c;
}

inline std::error_code make_error_code(Acts::EigenStepperError e) {
  return {static_cast<int>(e), Acts::EigenStepperErrorCategory()};
}
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::EigenStepperError> : std::true_type {};
}  // namespace std
