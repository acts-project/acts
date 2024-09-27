// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>

namespace Acts {

template <typename external_spacepoint_t>
class Seed {
 public:
  Seed(const external_spacepoint_t& b, const external_spacepoint_t& m,
       const external_spacepoint_t& u);

  void setVertexZ(float vertex);
  void setQuality(float seedQuality);

  const std::array<const external_spacepoint_t*, 3>& sp() const;
  float z() const;
  float seedQuality() const;

 private:
  std::array<const external_spacepoint_t*, 3> m_spacepoints;
  float m_vertexZ{0.f};
  float m_seedQuality{-std::numeric_limits<float>::infinity()};
};

}  // namespace Acts

#include "Acts/EventData/Seed.ipp"
