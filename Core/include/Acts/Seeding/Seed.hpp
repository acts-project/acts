// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

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
  inline bool operator<(const Seed& rhs) const {
    // not relevant ordering
    if (this->sp()[0]->z() != rhs.sp()[0]->z()) {
      return this->sp()[0]->z() < rhs.sp()[0]->z();
    } else if (this->sp()[1]->z() != rhs.sp()[1]->z()) {
      return this->sp()[1]->z() < rhs.sp()[1]->z();
    } else if (this->sp()[2]->z() != rhs.sp()[2]->z()) {
      return this->sp()[2]->z() < rhs.sp()[2]->z();
    }
    if (this->sp()[0]->x() != rhs.sp()[0]->x()) {
      return this->sp()[0]->x() < rhs.sp()[0]->x();
    } else if (this->sp()[1]->x() != rhs.sp()[1]->x()) {
      return this->sp()[1]->x() < rhs.sp()[1]->x();
    } else if (this->sp()[2]->x() != rhs.sp()[2]->x()) {
      return this->sp()[2]->x() < rhs.sp()[2]->x();
    }
    if (this->sp()[0]->y() != rhs.sp()[0]->y()) {
      return this->sp()[0]->y() < rhs.sp()[0]->y();
    } else if (this->sp()[1]->y() != rhs.sp()[1]->y()) {
      return this->sp()[1]->y() < rhs.sp()[1]->y();
    } else if (this->sp()[2]->y() != rhs.sp()[2]->y()) {
      return this->sp()[2]->y() < rhs.sp()[2]->y();
    }

    assert(this->sp()[0] == rhs.sp()[0]);
    assert(this->sp()[1] == rhs.sp()[1]);
    assert(this->sp()[2] == rhs.sp()[2]);
    return false;
  }

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
