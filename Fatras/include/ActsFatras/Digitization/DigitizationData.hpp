// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/EventData/ParameterSet.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>

#include <functional>

namespace Acts {
class Surface;
};  // namespace Acts

namespace ActsFatras {

class Hit;

/// Digitization input to be used by the digitizers to harmonize
/// the interface
///
struct DigitizationInput {
  std::reference_wrapper<const Hit> hit;
  std::reference_wrapper<const Acts::GeometryContext> geoContext;
  const Acts::Surface* surface = nullptr;

  /// Only valid constructor, wraps the @param hit_,
  /// the  and optionally the @param surface_
  DigitizationInput(
      std::reference_wrapper<const Hit> hit_,
      std::reference_wrapper<const Acts::GeometryContext> geoContext_,
      const Acts::Surface* surface_ = nullptr)
      : hit(hit_), geoContext(geoContext_), surface(surface_) {}

  DigitizationInput() = delete;
};

/// A single cell definition: index, value
using Cell = std::pair<unsigned int, double>;

/// A channel definition: Cell idenficiation, readout word
template <Acts::BoundIndices... kParameters>
using Channel = std::pair<std::array<Cell, sizeof...(kParameters)>, double>;

/// A Cluster definition: List of cells, parameter set, size
template <Acts::BoundIndices... kParameters>
struct Cluster {
  /// The parameters
  Acts::ParameterSet<Acts::BoundIndices, kParameters...> parameterSet;
  /// The resulting cluster size
  std::array<unsigned int, sizeof...(kParameters)> clusterSize;
  /// The contained Channels
  std::vector<const Channel<kParameters...> > channels;

  /// Cluster constructor
  Cluster(Acts::ParameterSet<Acts::BoundIndices, kParameters...> pSet,
          std::array<unsigned int, sizeof...(kParameters)> cSize,
          std::vector<const Channel<kParameters...> > cChannels)
      : parameterSet(pSet), clusterSize(cSize), channels(cChannels) {}
};

}  // namespace ActsFatras