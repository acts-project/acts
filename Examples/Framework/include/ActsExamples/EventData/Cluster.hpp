// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/Segmentizer.hpp"

#include <numeric>
#include <optional>
#include <vector>

namespace ActsExamples {

/// Simple struct holding cluster information.
struct Cluster {
  using Cell = ActsFatras::Segmentizer::ChannelSegment;
  std::size_t sizeLoc0 = 0;
  std::size_t sizeLoc1 = 0;
  std::vector<Cell> channels;

  // TODO make this be provided by Fatras?
  Acts::Vector3 globalPosition = Acts::Vector3::Zero();
  Acts::Vector3 localDirection = Acts::Vector3::Zero();
  Acts::Vector3 lengthDirection = Acts::Vector3::Zero();
  float localEta = 0.f;
  float localPhi = 0.f;
  float globalEta = 0.f;
  float globalPhi = 0.f;
  float etaAngle = 0.f;
  float phiAngle = 0.f;

  double sumActivations() const {
    return std::accumulate(
        channels.begin(), channels.end(), 0.0,
        [](double s, const Cluster::Cell& c) { return s + c.activation; });
  }
};

/// Clusters have a one-to-one relation with measurements
using ClusterContainer = std::vector<Cluster>;

}  // namespace ActsExamples
