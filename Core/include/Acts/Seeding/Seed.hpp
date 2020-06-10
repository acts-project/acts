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
template <typename spacepoint_t>
class Seed {
 public:
  Seed(const spacepoint_t& b, const spacepoint_t& m, const spacepoint_t& u,
       float vertex);
  Seed(const Seed&) = default;
  Seed& operator=(const Seed&);

  const std::vector<const spacepoint_t*>& sp() const { return m_spacepoints; }
  double z() const { return m_zvertex; }

 private:
  std::vector<const spacepoint_t*> m_spacepoints;
  float m_zvertex;
};

///////////////////////////////////////////////////////////////////////////////
// Constructors
///////////////////////////////////////////////////////////////////////////////

template <typename spacepoint_t>
Seed<spacepoint_t>::Seed(const spacepoint_t& b, const spacepoint_t& m,
                         const spacepoint_t& u, float vertex) {
  m_zvertex = vertex;
  m_spacepoints.push_back(&b);
  m_spacepoints.push_back(&m);
  m_spacepoints.push_back(&u);
}

}  // namespace Acts
