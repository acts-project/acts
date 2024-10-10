// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

struct SpacePoint {
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  int surface;
  float varianceR;
  float varianceZ;
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }
};

bool operator==(SpacePoint a, SpacePoint b) {
  if (fabs(a.m_x - b.m_x) < 1e-6 && fabs(a.m_y - b.m_y) < 1e-6 &&
      fabs(a.m_z - b.m_z) < 1e-6) {
    return true;
  } else {
    return false;
  }
}
