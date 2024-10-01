// This file is part of the ACTS project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Clusterization/Clusterization.hpp"
#include "Acts/Definitions/Algebra.hpp"

namespace Acts::Ccl {

template <typename Cell>
concept HasRetrievableTimeInfo = requires(Cell cell) {
  { getCellTime(cell) } -> std::same_as<Acts::ActsScalar>;
};

template <Acts::Ccl::HasRetrievableTimeInfo Cell, std::size_t N>
struct TimedConnect {};

}  // namespace Acts::Ccl
