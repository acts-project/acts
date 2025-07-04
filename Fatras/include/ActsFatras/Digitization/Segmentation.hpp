// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"

#include <array>

namespace ActsFatras {

/// Nested struct for stepping from one channel to the next.
struct ChannelStep {
  /// This is the delta to the last step in bins
  std::array<int, 2> delta = {0, 0};
  /// The intersection with the channel boundary
  Acts::Vector2 intersect;
  /// The patlength from the start
  double path = 0.;

  /// Constructor with arguments for a ChannelStep
  ///
  /// @param delta_ The bin delta for this step
  /// @param intersect_ The intersect with the channel boundary
  /// @param start The start of the surface segment, for path from origin
  ChannelStep(std::array<int, 2> delta_, Acts::Vector2 intersect_,
              const Acts::Vector2& start)
      : delta(delta_),
        intersect(std::move(intersect_)),
        path((intersect - start).norm()) {}

  /// Smaller operator for sorting the ChannelStep objects.
  ///
  /// @param cstep The other ChannelStep to be compared
  ///
  /// The ChannelStep objects can be compared with its path distance
  /// from the start (surface segment origin)
  bool operator<(const ChannelStep& cstep) const { return path < cstep.path; }
};

class ISegmentation {
 public:
  /// Virtual destructor
  virtual ~ISegmentation() = default;

  /// Return the bin according to the segmentation
  /// @param pos The position to be binned
  /// @return The bin in the segmentation
  virtual std::array<std::size_t, 2> bin(const Acts::Vector2& pos) const = 0;

  /// Return the position from the bin
  /// @param bin The bin to be converted
  ///
  /// @note no out of bounds checking done
  ///
  /// @return The position in the segmentation
  virtual Acts::Vector2 position(
      const std::array<std::size_t, 2>& bin) const = 0;

  /// The steps within the channels
  ///
  /// @param start The start of the segment in local 2D
  /// @param end The end of the segment in local 2D
  /// @return The channel steps
  virtual std::vector<ChannelStep> channelSteps(
      const Acts::Vector2& start, const Acts::Vector2& end) const = 0;
};

class CartesianSegmentation : public ISegmentation {
 public:
  CartesianSegmentation(const Acts::ProtoAxis& xAxis,
                        const Acts::ProtoAxis& yAxis);
  ~CartesianSegmentation() = default;

  /// @copydoc ISegmentation::bin
  std::array<std::size_t, 2> bin(const Acts::Vector2& pos) const override;

  /// @copydoc ISegmentation::position
  Acts::Vector2 position(const std::array<std::size_t, 2>& bin) const override;

  /// @copydoc ISegmentation::channelSteps
  std::vector<ChannelStep> channelSteps(
      const Acts::Vector2& start, const Acts::Vector2& end) const override;

 private:
  Acts::ProtoAxis m_xAxis;
  Acts::ProtoAxis m_yAxis;
};

}  // namespace ActsFatras
