// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>

namespace Acts {
namespace Experimental {

enum class Chi2FitterError {
  // ensure all values are non-zero
  NoMeasurementFound = 1,
};

std::error_code make_error_code(Chi2FitterError e);

}  // namespace Experimental
}  // namespace Acts

namespace std {
// register with STL
template <>
struct is_error_code_enum<Acts::Experimental::Chi2FitterError>
    : std::true_type {};
}  // namespace std
