// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>

namespace Acts {
namespace Experimental {

enum class GlobalChiSquareFitterError {
  // ensure all values are non-zero
  DetAIsZero = 1,
};

std::error_code make_error_code(
    Acts::Experimental::GlobalChiSquareFitterError e);

}  // namespace Experimental
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::Experimental::GlobalChiSquareFitterError>
    : std::true_type {};
}  // namespace std
