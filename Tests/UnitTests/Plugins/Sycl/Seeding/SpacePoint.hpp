// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

struct SpacePoint {
  SpacePoint()= default;;
  SpacePoint(float p_x, float p_y, float p_z, float p_r, int p_surface, float p_varianceR, float p_varianceZ):
    m_x(p_x),m_y(p_y),m_z(p_z),m_r(p_r),surface(p_surface),varianceR(p_varianceR),varianceZ(p_varianceZ){};
  float x() const { return m_x; }
  float y() const { return m_y; }
  float z() const { return m_z; }
  float r() const { return m_r; }
  friend bool operator==(const SpacePoint &a, const SpacePoint &b);
  float varianceR;
  float varianceZ;
  private:
    float m_x;
    float m_y;
    float m_z;
    float m_r;
    int surface;
};

bool operator==(const SpacePoint &a, const SpacePoint &b) {
  return  (a.m_x == b.m_x && a.m_y == b.m_y && a.m_z == b.m_z &&
      a.surface == b.surface && a.varianceR == b.varianceR &&
      a.varianceZ == b.varianceZ);
}