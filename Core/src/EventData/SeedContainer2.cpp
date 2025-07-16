// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SeedContainer2.hpp"

namespace Acts::Experimental {

void SeedContainer2::reserve(std::size_t size,
                             float averageSpacePoints) noexcept {
  m_sizes.reserve(size);
  m_spacePointOffsets.reserve(size);
  m_qualities.reserve(size);
  m_vertexZs.reserve(size);
  m_spacePoints.reserve(static_cast<std::size_t>(size * averageSpacePoints));
}

void SeedContainer2::clear() noexcept {
  m_sizes.clear();
  m_spacePointOffsets.clear();
  m_qualities.clear();
  m_vertexZs.clear();
  m_spacePoints.clear();
}

}  // namespace Acts::Experimental
