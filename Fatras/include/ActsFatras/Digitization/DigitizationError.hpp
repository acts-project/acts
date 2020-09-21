// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>        // for string printing
#include <system_error>  // bring in std::error_code et al

namespace ActsFatras {

enum class DigitizationError {
  SmearOutOfBounds = 1,
  SmearError = 2,
  NoSurfaceDefined = 3
};

namespace detail {
// /Define a custom error code category derived from std::error_category
class DigitizationErrorCategory : public std::error_category {
 public:
  /// Return a short descriptive name for the category
  const char* name() const noexcept final;

  /// Return what each enum means in text
  std::string message(int c) const final;
};
}  // namespace detail

/// Declare a global function returning a static instance of the custom category
const detail::DigitizationErrorCategory& DigitizationErrorCategory();

std::error_code make_error_code(DigitizationError e);

}  // namespace ActsFatras

namespace std {
/// register with STL
template <>
struct is_error_code_enum<ActsFatras::DigitizationError> : std::true_type {};
}  // namespace std
