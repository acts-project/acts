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
/// The current candidates are stored in a vector of candidates, where an index
/// is used to indicate the current active candidate.
///
/// @todo the NavigationStream should hold also the current volume it is in
/// if it represents the geometry stream.
struct NavigationStream {
  /// The query point for the navigation stream
  ///
  /// This holds the position and direction from which the navigation stream
  /// should either be initialized or updated.
  struct QueryPoint {
    /// The position of the query point
    Vector3 position = Vector3::Zero();
    /// The direction of the query point
    Vector3 direction = Vector3::Zero();
  };

  /// This is a candidate object of the navigation stream, it holds:
  ///
  /// - a Surface intersection
  /// - a Portal : set if the surface represents a portal
  /// - a BoundaryTolerance : the boundary tolerance used for the intersection
  struct Candidate {
    /// The intersection
    ObjectIntersection<Surface> intersection =
        ObjectIntersection<Surface>::invalid();
    /// The portal
    const Portal* portal = nullptr;
    /// The boundary tolerance
    BoundaryTolerance bTolerance = BoundaryTolerance::None();
    /// Target flag: true if this Candidate represents an abortTarget
    bool abortTarget = false;
    /// Convenience access to surface
    const Surface& surface() const { return *intersection.object(); }
    /// Cinvencience access to the path length
    ActsScalar pathLength() const { return intersection.pathLength(); }
  };

  /// The candidates of this navigation stream
  std::vector<Candidate> candidates;

  /// The currently active candidate
  std::size_t currentIndex = 0u;

  /// Switch to next next candidate
  ///
  /// @return true if a next candidate is available
  bool switchToNextCandidate() {
    if (currentIndex < candidates.size()) {
      ++currentIndex;
      return true;
    }
    return false;
  }

  /// Const access the current candidate
  const Candidate& currentCandidate() const {
    return candidates.at(currentIndex);
  }

  /// Non-cost access the current candidate
  Candidate& currentCandidate() { return candidates.at(currentIndex); }

  /// The number of active candidates
  std::size_t activeCandidates() const {
    return (candidates.size() - currentIndex);
  }
};

}  // namespace Acts
