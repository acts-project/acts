// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>

namespace Acts {

template <typename external_spacePoint_t>
class Seed {
 public:
  Seed(const external_spacePoint_t& b, const external_spacePoint_t& m,
       const external_spacePoint_t& u, float vertex,
       float seedQuality = -std::numeric_limits<float>::infinity());

  const std::array<const external_spacePoint_t*, 3>& sp() const;
  float z() const;
  float seedQuality() const;

 private:
  std::array<const external_spacePoint_t*, 3> m_spacepoints;
  float m_zvertex{0.};
  float m_seedQuality{-std::numeric_limits<float>::infinity()};
};

}  // namespace Acts

#include "Acts/EventData/Seed.ipp"
