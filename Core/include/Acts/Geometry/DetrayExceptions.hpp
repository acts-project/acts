// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

namespace Acts {
class DetrayNotAvailableException : public std::runtime_error {
 public:
  DetrayNotAvailableException();
};

inline DetrayNotAvailableException::DetrayNotAvailableException()
    : std::runtime_error("ACTS was built without detray support") {}

class DetrayUnsupportedMaterialException : public std::runtime_error {
 public:
  explicit DetrayUnsupportedMaterialException(std::string_view name);
};

inline DetrayUnsupportedMaterialException::DetrayUnsupportedMaterialException(
    std::string_view name)
    : std::runtime_error(std::string("Material type ") + std::string(name) +
                         " not supported by detray") {}

}  // namespace Acts
