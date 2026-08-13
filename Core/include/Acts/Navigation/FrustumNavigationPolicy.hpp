// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Frustum.hpp"

#include <numbers>

namespace Acts::Experimental {

/// A navigation policy that uses frustum-octree intersections to find portals

class FrustumNavigationPolicy : public INavigationPolicy {
 public:
  using BoundingBox = AxisAlignedBoundingBox<Volume, double, 3>;
  using Frustum3 = Frustum<double, 3, 3>;
  struct Config {
    /// The octree depth
    int depth = 3;
  };
  struct State {
    Frustum3 frustum =
        Frustum3(Vector3::Zero(), Vector3::Zero(), std::numbers::pi / 4);
  };

  /// Main constructor, which takes the top-level volume and builds the octree
  /// @param gctx The geometrycontext object
  /// @param volume The tracking volume used to build the octree
  /// @param config The configuration of the Navigation Policy
  /// @param logger A logging instance
  explicit FrustumNavigationPolicy(const GeometryContext& gctx,
                                   const TrackingVolume& volume,
                                   const Logger& logger, const Config& config);

  /// Update the navigation state
  /// @param gctx The geometry context
  /// @param args The navigation arguments
  /// @param state The navigation state for this policy
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const GeometryContext& gctx,
                            const NavigationArguments& /*args*/,
                            NavigationPolicyState& state,
                            AppendOnlyNavigationStream& stream,
                            const Logger& logger) const;

  /// Connect this policy with a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override;

  /// Check the validity of the navigation state
  /// @param args The navigation arguments
  /// @param state The navigation state to check
  /// @param logger The logger
  bool isValid(const GeometryContext& /*gctx*/, const NavigationArguments& args,
               NavigationPolicyState& state,
               const Logger& logger) const override;

  /// Create and initialize a new navigation state
  /// @param args The navigation arguments
  /// @param stateManager The navigation state manager to push the new state onto
  /// @param logger The logger
  void createState(const GeometryContext& /*gctx*/,
                   const NavigationArguments& args,
                   NavigationPolicyStateManager& stateManager,
                   const Logger& logger) const override;

  /// Remove the state from the state manager
  /// @param stateManager The state manager to pop the state from
  /// @param logger The logger
  void popState(NavigationPolicyStateManager& stateManager,
                const Logger& logger) const override;

 private:
  // The vector of unique pointers of bounding boxes
  std::vector<std::unique_ptr<BoundingBox> > m_boxes;

  // The top-level bounding box with the octree
  BoundingBox* m_topBox;

  // associated volume name
  std::string m_name;
};

static_assert(NavigationPolicyConcept<FrustumNavigationPolicy>);

}  // namespace Acts::Experimental
