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
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Seeding/GbtsGeometry.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>

namespace Acts::Experimental {

struct GbtsEdgeState {

public:

struct Compare {
    bool operator()(const struct GbtsEdgeState* s1, const struct GbtsEdgeState* s2) {
      return s1->m_J > s2->m_J;
    }
  };


  GbtsEdgeState() {};

  GbtsEdgeState(bool f) : m_initialized(f) {};

  ~GbtsEdgeState() {};

  void initialize(GbtsEdge*);
  void clone(const struct GbtsEdgeState&);

  float m_J{};

  std::vector<GbtsEdge*> m_vs;

  float m_X[3]{}, m_Y[2]{}, m_Cx[3][3]{}, m_Cy[2][2]{};
  float m_refX{}, m_refY{}, m_c{}, m_s{};
  
  bool m_initialized{false};

};

#define MAX_EDGE_STATE 2500

class GbtsTrackingFilter {
 public:
  GbtsTrackingFilter(const std::vector<TrigInDetSiLayer>&, std::vector<GbtsEdge>&, SeedFinderGbtsConfig& config);
  ~GbtsTrackingFilter(){};

  void followTrack(GbtsEdge*, GbtsEdgeState&);

 protected:

  void propagate(GbtsEdge*, GbtsEdgeState&);

  bool update(GbtsEdge*, GbtsEdgeState&);

  int getLayerType(int);  


  const std::vector<TrigInDetSiLayer>& m_geo;
  
  std::vector<GbtsEdge>& m_segStore;
 
  std::vector<GbtsEdgeState*> m_stateVec;

  GbtsEdgeState m_stateStore[MAX_EDGE_STATE];

  int m_globalStateCounter{0};

  SeedFinderGbtsConfig& m_config;

};



}  // namespace Acts::Experimental
