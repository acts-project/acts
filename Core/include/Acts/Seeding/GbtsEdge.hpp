// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/Types.hpp"

#include <array>

namespace Acts::Experimental {

//! [gbts edge]
/// Edge between two GBTS nodes with fit parameters.
struct GbtsEdge final {
  GbtsEdge() = default;

  /// Constructor
  /// @param n1_ Inner node index
  /// @param n2_ Outer node index
  /// @param tauRatioCut_ Tau ratio cut of this doublet as an outer pair
  /// @param n2PixelBarrel_ Whether the outer node is on the pixel barrel
  /// @param p1_ First fit parameter
  /// @param p2_ Second fit parameter
  /// @param p3_ Third fit parameter
  GbtsEdge(SpacePointIndex n1_, SpacePointIndex n2_, float tauRatioCut_,
           bool n2PixelBarrel_, float p1_, float p2_, float p3_)
      : n1{n1_},
        n2{n2_},
        n2PixelBarrel{n2PixelBarrel_},
        p{p1_, p2_, p3_},
        tauRatioCut{tauRatioCut_} {}

  /// Inner node of the edge
  SpacePointIndex n1{kSpacePointIndexInvalid};
  /// Outer node of the edge
  SpacePointIndex n2{kSpacePointIndexInvalid};

  /// Whether the outer node is on the pixel barrel.
  bool n2PixelBarrel{};

  std::array<float, 3> p{};

  /// GbtsLayerConnection::tauRatioCutOuter of this doublet's layer pair,
  /// cached for the neighbour loop.
  float tauRatioCut{};
};
//! [gbts edge]

}  // namespace Acts::Experimental
