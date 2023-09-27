// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>

namespace Acts {
namespace Experimental {

enum class GsfError {
  NoMeasurementStatesCreatedForward = 1,
  NoMeasurementStatesCreatedBackward,
  NoMeasurementStatesCreatedFinal,
  StartParametersNotOnStartSurface
};

std::error_code make_error_code(GsfError e);

}  // namespace Experimental
}  // namespace Acts

// register with STL
namespace std {
template <>
struct is_error_code_enum<Acts::Experimental::GsfError> : std::true_type {};
}  // namespace std
