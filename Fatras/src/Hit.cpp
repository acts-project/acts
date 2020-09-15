// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsFatras/EventData/Hit.hpp"

ActsFatras::Hit::Hit(Acts::GeometryIdentifier geometryId, Barcode particleId,
                     const Vector4& pos4, const Vector4& before4,
                     const Vector4& after4, int32_t index_)
    : m_geometryId(geometryId),
      m_particleId(particleId),
      m_index(index_),
      m_pos4(pos4),
      m_before4(before4),
      m_after4(after4) {}

bool ActsFatras::Hit::operator==(const ActsFatras::Hit& other) const {
  if (this != &other) {
    bool identical =
        (m_geometryId == other.m_geometryId and
         m_particleId == other.m_particleId and m_index == other.m_index);
    if (identical) {
      return (m_pos4 == other.m_pos4 and m_before4 == other.m_before4 and
              m_after4 == other.m_after4);
    }
    return false;
  }
  return true;
}