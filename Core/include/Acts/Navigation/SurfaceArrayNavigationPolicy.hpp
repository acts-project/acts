// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Navigation/INavigationPolicy.hpp"

#pragma once

namespace Acts {

class SurfaceArray;

/// A navigation policy that internally uses the Gen1 @c SurfaceArray class
class SurfaceArrayNavigationPolicy : public INavigationPolicy {
 public:
  /// Enum for configuring which type of surface array to produce. This affects
  /// the projection that is used for creating the binning structure.
  enum class LayerType { Cylinder, Disc, Plane };

  /// Config struct to configure the surface array navigation
  struct Config {
    /// The type of the layer
    LayerType layerType = LayerType::Cylinder;
    /// The number of bins in the local directions. The interpretation depends
    /// on the layer type.
    std::pair<std::size_t, std::size_t> bins;
  };

  /// Main constructor, which internally creates the surface array acceleration
  /// structure.
  /// @note Expects that all relevant surfaces are registered with @p volume.
  ///       Only selects sensitive surfaces for the surface array.
  /// @param gctx The geometry context
  /// @param volume The *layer volume* to construct the surface array from
  /// @param logger A logging instance
  /// @param config The configuration for the surface array
  explicit SurfaceArrayNavigationPolicy(const GeometryContext& gctx,
                                        const TrackingVolume& volume,
                                        const Logger& logger, Config config);

  /// Update the navigation state from the surface array
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation policy state
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& args,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect this policy with a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override;

  /// Output stream operator for the contained layer type enum
  /// @param os The output stream
  /// @param layerType The layer type to print
  friend std::ostream& operator<<(std::ostream& os,
                                  const LayerType& layerType) {
    switch (layerType) {
      using enum LayerType;
      case Cylinder:
        os << "Cylinder";
        break;
      case Disc:
        os << "Disc";
        break;
      case Plane:
        os << "Plane";
        break;
    }
    return os;
  }

  /// Const reference access to the layer array
  /// @return The surface array
  const SurfaceArray& surfaceArray() const;

 private:
  std::unique_ptr<SurfaceArray> m_surfaceArray{};
  const TrackingVolume& m_volume;
};

static_assert(NavigationPolicyConcept<SurfaceArrayNavigationPolicy>);

}  // namespace Acts
