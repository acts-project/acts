// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SeedContainer2.hpp"

namespace Acts::Experimental {

void SeedContainer2::reserve(std::size_t size, float averageSpacePoints) {
  m_entries.reserve(size);
  m_spacePoints.reserve(static_cast<std::size_t>(size * averageSpacePoints));
}

void SeedContainer2::clear() {
  m_entries.clear();
  m_spacePoints.clear();
}

}  // namespace Acts::Experimental
