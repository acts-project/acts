// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts::detail {

/// This helper function calculates the combination
/// of two phi sectors, defined by a phi half-length +
/// a half phi sector in the range [0,pi). The two
/// ranges need to line up, i.e. that one of the sector
/// ends exactly where the other one starts.
std::tuple<double, double, bool> mergedPhiSector(
    double hlPhi1, double avgPhi1, double hlPhi2, double avgPhi2,
    const Logger& logger = getDummyLogger(),
    double tolerance = s_onSurfaceTolerance);

}  // namespace Acts::detail
