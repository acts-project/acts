// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <memory>

namespace Acts {
template <typename SpacePoint>
class InternalSeed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  InternalSeed(InternalSpacePoint<SpacePoint>& s0,
               InternalSpacePoint<SpacePoint>& s1,
               InternalSpacePoint<SpacePoint>& s2, float z,
               bool qualitySeed = false);
  InternalSeed& operator=(const InternalSeed& seed);

  const std::array<InternalSpacePoint<SpacePoint>*, 3> sp;
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

template <typename SpacePoint>
inline InternalSeed<SpacePoint>& InternalSeed<SpacePoint>::operator=(
    const InternalSeed<SpacePoint>& seed) {
  m_z = seed.m_z;
  sp = seed.sp;
  m_qualitySeed = seed.m_qualitySeed;
  return (*this);
}

template <typename SpacePoint>
inline InternalSeed<SpacePoint>::InternalSeed(
    InternalSpacePoint<SpacePoint>& s0, InternalSpacePoint<SpacePoint>& s1,
    InternalSpacePoint<SpacePoint>& s2, float z, bool qualitySeed)
    : sp({&s0, &s1, &s2}) {
  m_z = z;
  m_qualitySeed = qualitySeed;
}

/// @endcond

}  // namespace Acts
