// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

struct SpacePoint {
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  int surface;
  float varianceR;
  float varianceZ;
  float m_topHalfStripLength = 0;
  float m_bottomHalfStripLength = 0;
  Acts::Vector3 m_topStripDirection = {0, 0, 0};
  Acts::Vector3 m_bottomStripDirection = {0, 0, 0};
  Acts::Vector3 m_stripCenterDistance = {0, 0, 0};
  Acts::Vector3 m_bottomStripCenterPosition = {0, 0, 0};
  bool m_validDoubleMeasurementDetails = false;
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }
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
};

bool operator==(SpacePoint a, SpacePoint b) {
  if (fabs(a.m_x - b.m_x) < 1e-6 && fabs(a.m_y - b.m_y) < 1e-6 &&
      fabs(a.m_z - b.m_z) < 1e-6) {
    return true;
  } else {
    return false;
  }
}
