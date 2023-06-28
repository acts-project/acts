// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename external_spacePoint_t>
Seed<external_spacePoint_t>::Seed(const external_spacePoint_t& b,
                                  const external_spacePoint_t& m,
                                  const external_spacePoint_t& u, float vertex,
                                  float seedQuality)
    : m_spacepoints({&b, &m, &u}),
      m_zvertex(vertex),
      m_seedQuality(seedQuality) {}

template <typename external_spacePoint_t>
inline const std::array<const external_spacePoint_t*, 3>&
Seed<external_spacePoint_t>::sp() const {
  return m_spacepoints;
}

template <typename external_spacePoint_t>
inline float Seed<external_spacePoint_t>::z() const {
  return m_zvertex;
}

template <typename external_spacePoint_t>
inline float Seed<external_spacePoint_t>::seedQuality() const {
  return m_seedQuality;
}

}  // namespace Acts
