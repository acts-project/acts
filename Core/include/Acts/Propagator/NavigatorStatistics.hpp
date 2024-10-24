// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>

namespace Acts {

struct NavigatorStatistics {
  std::size_t nRenavigations = 0;

  std::size_t nVolumeSwitches = 0;

  std::size_t nBoundaryCandidates = 0;
  std::size_t nBoundaryHits = 0;
  std::size_t nBoundaryDiscards = 0;

  std::size_t nLayerCandidates = 0;
  std::size_t nLayerHits = 0;
  std::size_t nLayerDiscards = 0;

  std::size_t nSurfaceCandidates = 0;
  std::size_t nSurfaceHits = 0;
  std::size_t nSurfaceDiscards = 0;

  std::size_t nTotalCandidates = 0;
  std::size_t nTotalHits = 0;
  std::size_t nTotalDiscards = 0;
};

}  // namespace Acts
