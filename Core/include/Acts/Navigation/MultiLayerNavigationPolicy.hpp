// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"


#pragma once

namespace Acts {

/// A navigation policy that uses SurfaceArray to create surfaces on a plane layer 
/// Navigate through a multilayer structure by creating an artificial on the grid.
class MultiLayerNavigationPolicy : public SurfaceArrayNavigationPolicy {
  public:
  ///Constructor from the SurfaceArrayNavigation
  using SurfaceArrayNavigationPolicy::SurfaceArrayNavigationPolicy;


  /// Update the navigation state from the surface array
  /// @param args The navigation arguments
  /// @param stream The navigation stream to update
  /// @param logger The logger
  void initializeCandidates(const NavigationArguments& args,
    AppendOnlyNavigationStream& stream,
    const Logger& logger) const override final;

  /// Connect this policy with a navigation delegate
  /// @param delegate The navigation delegate to connect to
  void connect(NavigationDelegate& delegate) const override;

  private:
  /// Generate a path in the multilayer
  /// @param startPosition The starting position of the path (in local frame)
  /// @param direction The direction of the path (in local frame)
  /// @param stepSize The step size for the path
  /// @param numberOfSteps The number of steps to take
  /// @return A vector of positions along the path
  std::vector<Vector3> generatePath(const Vector3& startPosition,
    const Vector3& direction, double stepSize,
    std::size_t numberOfSteps) const;


};

static_assert(NavigationPolicyConcept<MultiLayerNavigationPolicy>);

} // namespace Acts
