// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style
#include "Acts/Seeding/GbtsDataStorage.hpp"  //includes geo which has trigindetsilayer, may move this to trigbase
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"

#include <cmath>
#include <cstring>
#include <vector>

namespace Acts::Experimental {

struct GbtsEdgeState {
 public:
  struct Compare {
    bool operator()(const struct GbtsEdgeState* s1,
                    const struct GbtsEdgeState* s2) {
      return s1->m_J > s2->m_J;
    }
  };

  GbtsEdgeState() = default;

  explicit GbtsEdgeState(bool f) : m_initialized(f) {}

  void initialize(GbtsEdge& pS);
  void clone(const struct GbtsEdgeState& st);

  float m_J{};

  std::vector<GbtsEdge*> m_vs;

  float m_X[3]{}, m_Y[2]{}, m_Cx[3][3]{}, m_Cy[2][2]{};
  float m_refX{}, m_refY{}, m_c{}, m_s{};

  bool m_initialized{false};
};

class GbtsTrackingFilter {
 public:
  static constexpr std::int32_t MAX_EDGE_STATE = 2500;
  GbtsTrackingFilter(const std::vector<TrigInDetSiLayer>& g,
                     std::vector<GbtsEdge>& sb,
                     const SeedFinderGbtsConfig& config);

  void followTrack(GbtsEdge& pS, GbtsEdgeState& output);

 protected:
  void propagate(GbtsEdge& pS, GbtsEdgeState& ts);

  bool update(GbtsEdge& pS, GbtsEdgeState& ts);

  std::int32_t getLayerType(std::int32_t l);

  const std::vector<TrigInDetSiLayer>& m_geo;

  std::vector<GbtsEdge>& m_segStore;

  std::vector<GbtsEdgeState*> m_stateVec;

  std::array<GbtsEdgeState, MAX_EDGE_STATE> m_stateStore;

  std::int32_t m_globalStateCounter{0};

  const SeedFinderGbtsConfig& m_config;
};

}  // namespace Acts::Experimental
