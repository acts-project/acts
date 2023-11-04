// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GNN_TrackingFilter.hpp"

#include <cmath>

#define MAX_SILICON_LAYER_NUM 19
#define OffsetEndcapPixels 7
#define OffsetBarrelSCT 3
#define OffsetEndcapSCT 10

template <typename space_point_t>
class TrigInDetTriplet {
 public:
  TrigInDetTriplet() = delete;  // to prevent creation w/o initialization

  TrigInDetTriplet(Acts::FTF_SP<space_point_t> s1,
                   Acts::FTF_SP<space_point_t> s2,
                   Acts::FTF_SP<space_point_t> s3, float Q)
      : m_s1(std::move(s1)), m_s2(std::move(s2)), m_s3(std::move(s3)), m_Q(Q) {}

  TrigInDetTriplet(TrigInDetTriplet* t)
      : m_s1(t->m_s1), m_s2(t->m_s2), m_s3(t->m_s3), m_Q(t->m_Q) {}

  const Acts::FTF_SP<space_point_t>& s1() const { return m_s1; }
  const Acts::FTF_SP<space_point_t>& s2() const { return m_s2; }
  const Acts::FTF_SP<space_point_t>& s3() const { return m_s3; }
  float Q() const { return m_Q; }
  void Q(double newQ) { m_Q = newQ; }

 protected:
  Acts::FTF_SP<space_point_t> m_s1;
  Acts::FTF_SP<space_point_t> m_s2;
  Acts::FTF_SP<space_point_t> m_s3;
  float m_Q;  // Quality
};
