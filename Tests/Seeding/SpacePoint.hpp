// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
struct SpacePoint{
  float m_x;
  float m_y;
  float m_z;
  float m_r;
  int surface;
  float covr;
  float covz;
  float x() const {
    return m_x;
  }
  float y() const {
    return m_y;
  }
  float z() const {
    return m_z;
  }
  float r() const {
    return m_r;
  }
};
