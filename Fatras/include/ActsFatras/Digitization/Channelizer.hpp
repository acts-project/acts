// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/Digitization/DigitizationData.hpp"
#include "ActsFatras/Digitization/detail/WeightedChannelCombiner.hpp"
#include <Acts/EventData/ParameterSet.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/SurfaceError.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <Acts/Utilities/ParameterDefinitions.hpp>
#include <array>

namespace ActsFatras {

/// The Channelizer splits a full surface segment (i.e. after drift/projects)
/// into segements
///
class Channelizer {
 public:
  /// Nexted struct for stepping from one channel to the next.
  struct ChannelStep {
    std::array<int, 2> delta = {0, 0};
    Acts::Vector2D intersect;
    double path = 0.;

    /// Constructor with arguments
    ///
    /// @param delta_ The bin delta step
    /// @param intersect_ The intersect with the channel boundary
    /// @param start The start of the surface segment
    ChannelStep(std::array<int, 2> delta_, Acts::Vector2D intersect_,
                const Acts::Vector2D& start)
        : delta(delta_),
          intersect(intersect_),
          path((intersect - start).norm()) {}

    /// Smaller operator for sorting the ChannelSteps
    bool operator<(const ChannelStep& cstep) const { return path < cstep.path; }
  };

  /// Nested struct for representing channel steps.
  struct ChannelSegment {
    std::array<unsigned int, 2> bin = {0, 0};
    double value = 0.;

    /// Constructor with arguments
    ///
    /// @param bin_ The bin corresponding to this step
    /// @param value_ The segment length for this bin
    ChannelSegment(std::array<unsigned int, 2> bin_, double value_)
        : bin(std::move(bin_)), value(value_) {}
  };

  /// Divide the surface segment into channel segments.
  ///
  /// @param dInput The digitization input including segmentation
  /// @param start The surface segment start (cartesian coordinates)
  /// @param send The surface segement end (cartesian coordinates)
  std::vector<ChannelSegment> segments(const DigitizationInput& dInput,
                                       const Acts::Vector2D& start,
                                       const Acts::Vector2D& end) const;
};

}  // namespace ActsFatras
