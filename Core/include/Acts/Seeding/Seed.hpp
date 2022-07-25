// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SpacePoint.hpp"

#include <memory>

namespace Acts {
class InternalSeed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSeed(Acts::SpacePoint& s0, Acts::SpacePoint& s1, Acts::SpacePoint& s2,
               float z, bool qualitySeed = false);

  std::array<Acts::SpacePoint*, 3> sp;
  float z() const { return m_z; }
  bool qualitySeed() const { return m_qualitySeed; }

 protected:
  float m_z;
  bool m_qualitySeed;
};

/// @cond

/////////////////////////////////////////////////////////////////////////////////
// Inline methods
/////////////////////////////////////////////////////////////////////////////////

inline InternalSeed::InternalSeed(Acts::SpacePoint& s0, Acts::SpacePoint& s1,
                                  Acts::SpacePoint& s2, float z,
                                  bool qualitySeed)
    : sp({&s0, &s1, &s2}) {
  m_z = z;
  m_qualitySeed = qualitySeed;
}

using SeedContainer = std::vector<Acts::InternalSeed>;
/// @endcond

}  // namespace Acts
