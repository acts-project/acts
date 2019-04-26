// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

namespace Acts {
template <typename SpacePoint>
class Seed {
  /////////////////////////////////////////////////////////////////////////////////
  // Public methods:
  /////////////////////////////////////////////////////////////////////////////////

 public:
  Seed(const SpacePoint& b, const SpacePoint& m, const SpacePoint& u,
       float vertex);
  Seed(const Seed&) = default;
  Seed& operator=(const Seed&);

  const std::vector<const SpacePoint*>& sp() const { return m_spacepoints; }
  double z() const { return m_zvertex; }

 private:
  std::vector<const SpacePoint*> m_spacepoints;
  float m_zvertex;
};

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

template <typename SpacePoint>
Seed<SpacePoint>::Seed(const SpacePoint& b, const SpacePoint& m,
                       const SpacePoint& u, float vertex) {
  m_zvertex = vertex;
  m_spacepoints.push_back(&b);
  m_spacepoints.push_back(&m);
  m_spacepoints.push_back(&u);
}

}  // namespace Acts
