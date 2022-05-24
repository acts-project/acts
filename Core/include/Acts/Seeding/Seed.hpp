// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <limits>
#include <vector>

#include <boost/container/small_vector.hpp>

namespace Acts {
template <typename SpacePoint>
class Seed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  Seed(const SpacePoint& b, const SpacePoint& m, const SpacePoint& u,
       float vertex,
       float seedQuality = -std::numeric_limits<float>::infinity());
  Seed(const Seed&) = default;
  Seed& operator=(const Seed&) = default;

  const auto& sp() const { return m_spacepoints; }
  double z() const { return m_zvertex; }
  float seedQuality() const { return m_seedQuality; }

 private:
  boost::container::small_vector<const SpacePoint*, 3> m_spacepoints;
  float m_zvertex;
  float m_seedQuality;
};

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
Seed<SpacePoint>::Seed(const SpacePoint& b, const SpacePoint& m,
                       const SpacePoint& u, float vertex, float seedQuality) {
  m_zvertex = vertex;
  m_spacepoints.push_back(&b);
  m_spacepoints.push_back(&m);
  m_spacepoints.push_back(&u);
  m_seedQuality = seedQuality;
}

}  // namespace Acts
