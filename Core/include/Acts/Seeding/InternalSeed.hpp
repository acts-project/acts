// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <memory>

namespace Acts {

/// @brief The Internal seed representation, consisntent of 3 space points
/// and a longitudinal parameter
template <typename external_spacepoint_t>
class InternalSeed {
 public:
  /// Constructor with arguments
  /// @param s0 first space point
  /// @param s1 second space point
  /// @param s2 third space point
  /// @param z longitudinal impact
  InternalSeed(const InternalSpacePoint<external_spacepoint_t>& s0,
               const InternalSpacePoint<external_spacepoint_t>& s1,
               const InternalSpacePoint<external_spacepoint_t>& s2, float z);

  InternalSeed& operator=(const InternalSeed& seed);

  const std::array<const InternalSpacePoint<external_spacepoint_t>*, 3> sp;
  float z() const { return m_z; }

 protected:
  float m_z;
};

template <typename external_spacepoint_t>
inline InternalSeed<external_spacepoint_t>&
InternalSeed<external_spacepoint_t>::operator=(
    const InternalSeed<external_spacepoint_t>& seed) {
  m_z = seed.m_z;
  sp = seed.sp;
  return (*this);
}

template <typename external_spacepoint_t>
inline InternalSeed<external_spacepoint_t>::InternalSeed(
    const InternalSpacePoint<external_spacepoint_t>& s0,
    const InternalSpacePoint<external_spacepoint_t>& s1,
    const InternalSpacePoint<external_spacepoint_t>& s2, float z)
    : sp({&s0, &s1, &s2}) {
  m_z = z;
}

}  // namespace Acts
