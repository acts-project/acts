// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

struct SpacePoint {
  SpacePoint() = default;
  ;
  SpacePoint(float p_x, float p_y, float p_z, float p_r, int p_surface,
             float p_varianceR, float p_varianceZ, int p_id)
      : varianceR(p_varianceR),
        varianceZ(p_varianceZ),
        m_surface(p_surface),
        m_id(p_id),
        m_x(p_x),
        m_y(p_y),
        m_z(p_z),
        m_r(p_r){};

  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }
  int id() const { return m_id; }
  int surface() const { return m_surface; }
  constexpr float topHalfStripLength() const { return m_topHalfStripLength; }
  constexpr float bottomHalfStripLength() const {
    return m_bottomHalfStripLength;
  }
  Acts::Vector3 topStripDirection() const { return m_topStripDirection; }
  Acts::Vector3 bottomStripDirection() const { return m_bottomStripDirection; }
  Acts::Vector3 stripCenterDistance() const { return m_stripCenterDistance; }
  Acts::Vector3 bottomStripCenterPosition() const {
    return m_bottomStripCenterPosition;
  }
  constexpr bool validDoubleMeasurementDetails() const {
    return m_validDoubleMeasurementDetails;
  }
  friend bool operator==(const SpacePoint &a, const SpacePoint &b);
  float varianceR;
  float varianceZ;
  int m_surface;

 private:
  int m_id;
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  float m_topHalfStripLength = 0;
  float m_bottomHalfStripLength = 0;
  Acts::Vector3 m_topStripDirection = {0, 0, 0};
  Acts::Vector3 m_bottomStripDirection = {0, 0, 0};
  Acts::Vector3 m_stripCenterDistance = {0, 0, 0};
  Acts::Vector3 m_bottomStripCenterPosition = {0, 0, 0};
  bool m_validDoubleMeasurementDetails = false;
};

bool operator==(const SpacePoint &a, const SpacePoint &b) {
  return a.m_id == b.m_id;
}
