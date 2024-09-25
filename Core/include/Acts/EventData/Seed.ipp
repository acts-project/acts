// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename external_spacepoint_t>
Seed<external_spacepoint_t>::Seed(const external_spacepoint_t& b,
                                  const external_spacepoint_t& m,
                                  const external_spacepoint_t& u)
    : m_spacepoints({&b, &m, &u}) {}

template <typename external_spacepoint_t>
void Seed<external_spacepoint_t>::setVertexZ(float vertex) {
  m_vertexZ = vertex;
}

template <typename external_spacepoint_t>
void Seed<external_spacepoint_t>::setQuality(float seedQuality) {
  m_seedQuality = seedQuality;
}

template <typename external_spacepoint_t>
const std::array<const external_spacepoint_t*, 3>&
Seed<external_spacepoint_t>::sp() const {
  return m_spacepoints;
}

template <typename external_spacepoint_t>
float Seed<external_spacepoint_t>::z() const {
  return m_vertexZ;
}

template <typename external_spacepoint_t>
float Seed<external_spacepoint_t>::seedQuality() const {
  return m_seedQuality;
}

}  // namespace Acts
