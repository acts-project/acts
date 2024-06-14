// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <utility>

namespace Acts::detail {

/// This helper functions calculates the combination
/// of two phi sectors, defined by a phi half-length +
/// a half phi sector in the range [0,pi). The two
/// ranges need to line up, i.e. that one of the sector
/// ends exactly where the other one starts.
std::pair<ActsScalar, ActsScalar> mergedPhiSector(
    ActsScalar hlPhi1, ActsScalar avgPhi1, ActsScalar hlPhi2,
    ActsScalar avgPhi2, const Logger& logger = getDummyLogger(),
    ActsScalar tolerance = s_onSurfaceTolerance);

}  // namespace Acts::detail
