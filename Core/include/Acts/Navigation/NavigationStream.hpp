// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <limits>
#include <tuple>
#include <vector>

namespace Acts {

// To be removed when the namespace Experimental is omitted
namespace Experimental {
class Portal;
}
using namespace Experimental;

class Surface;

/// @brief The NavigationStream is a container for the navigation candidates
///
/// The current candidates are stored in a vector with a range defined
/// by a pair of indices. This implementation allows passed or unreachable
/// candidates to be shaddowed without removing them from the container.
///
/// @todo the NavigationStream should hold also the current volume it is in
/// if it represents the geometry stream.
///
/// A surface proximity parameter can be used to chose which sort of
/// intersection path length update is needed.
struct NavigationStream {
  /// @brief the Query ppoint for updating the navigation stream
  struct QueryPoint {
    /// The position of the query point
    Vector3 position = Vector3::Zero();
    /// The direction of the query point
    Vector3 direction = Vector3::Zero();
  };

  /// This is a candidate type for a Surface intersection
  using SurfaceIntersection = ObjectIntersection<Surface>;

  /// This is a candidate object of the navigation stream
  /// a Surface intersection
  /// a Portal : set if the surface represents a portal
  /// a BoundaryTolerance : the boundary tolerance used for the intersection
  struct Candidate {
    /// The intersection
    SurfaceIntersection intersection = SurfaceIntersection::invalid();
    /// The portal
    const Portal* portal = nullptr;
    /// The boundary tolerance
    BoundaryTolerance bTolerance = BoundaryTolerance::None();
    /// Convenience access to surface
    const Surface& surface() const { return *intersection.object(); }
    /// Cinvencience access to the path length
    ActsScalar pathLength() const { return intersection.pathLength(); }
  };

  /// The candidates for the navigation
  std::vector<Candidate> candidates;

  /// The currently active candidate range
  size_t currentIndex = 0u;

  /// Progress to next next candidate
  ///
  /// @return true if a next candidate is available
  bool switchToNextCandidate() {
    if (currentIndex < candidates.size()) {
      ++currentIndex;
      return true;
    }
    return false;
  }

  /// @brief Const-access the current candidate
  const Candidate& currentCandidate() const {
    return candidates.at(currentIndex);
  }

  /// @brief Noncost-access the current candidate
  Candidate& currentCandidate() { return candidates.at(currentIndex); }

  /// @brief Access the current candidate
  std::size_t activeCandidates() const {
    return (candidates.size() - currentIndex);
  }
};

}  // namespace Acts
