// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/EventData/ParameterSet.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Utilities/BinUtility.hpp>

#include <functional>
#include <unordered_set>

namespace Acts {
class Surface;
}  // namespace Acts

namespace ActsFatras {

class Hit;

/// Digitization input to be used by the digitizers to harmonize
/// the interface
///
struct DigitizationInput {
  std::reference_wrapper<const Hit> hit;
  std::reference_wrapper<const Acts::GeometryContext> geoContext;
  const Acts::Surface* surface = nullptr;
  Acts::BinUtility segmentation;

  /// Only valid constructor, wraps the @param hit_,
  /// the  and optionally the @param surface_
  DigitizationInput(
      std::reference_wrapper<const Hit> hit_,
      std::reference_wrapper<const Acts::GeometryContext> geoContext_,
      const Acts::Surface* surface_ = nullptr,
      Acts::BinUtility segmentation_ = Acts::BinUtility())
      : hit(hit_),
        geoContext(geoContext_),
        surface(surface_),
        segmentation(segmentation_) {}

  DigitizationInput() = delete;
};

/// A single cell definition: index, cell central value
using Cell = std::pair<unsigned int, double>;

/// A channel definition: Cell identification, readout word, links
///
/// @tparam signal_t Type of the signal, requires += operator
/// @tparam kParameters ... The parameters pack
template <typename signal_t, Acts::BoundIndices... kParameters>
struct Channel {
  /// The cell identification in sizeof..(kParameters) dimensions
  std::array<Cell, sizeof...(kParameters)> cellId;
  /// The signal value, as complex as possible,
  /// but need += operator and double() cast for the weight
  signal_t value = 0.;
  /// The potential (truth) links
  std::unordered_set<unsigned int> links = {};

  /// Channel constructor
  ///
  /// @param cellId_ The Cell idenficiation and position
  /// @param value_ The Cell value
  /// @param links_ The (optional) links to e.g. truth indices
  Channel(std::array<Cell, sizeof...(kParameters)> cellId_, signal_t value_,
          std::unordered_set<unsigned int> links_ = {})
      : cellId(cellId_), value(value_), links(links_) {}

  Channel() = delete;
};

/// A Cluster definition.
///
/// @tparam signal_t Type of the signal carried, see above
/// @tparam kParameters Parameters pack for the cluster
///
template <typename signal_t, Acts::BoundIndices... kParameters>
struct Cluster {
  /// The parameters
  Acts::ParameterSet<Acts::BoundIndices, kParameters...> parameterSet;
  /// The resulting cluster size
  std::array<unsigned int, sizeof...(kParameters)> clusterSize;
  /// The contained Channels
  std::vector<Channel<signal_t, kParameters...> > channels;

  /// Cluster constructor
  ///
  /// @param pSet The parameter set repesenting cluster position and error
  /// @param cSize The cluster size definition
  /// @param cChannels The channel
  Cluster(Acts::ParameterSet<Acts::BoundIndices, kParameters...> pSet,
          std::array<unsigned int, sizeof...(kParameters)> cSize,
          std::vector<Channel<signal_t, kParameters...> > cChannels)
      : parameterSet(pSet), clusterSize(cSize), channels(cChannels) {}

  Cluster() = delete;
};

}  // namespace ActsFatras