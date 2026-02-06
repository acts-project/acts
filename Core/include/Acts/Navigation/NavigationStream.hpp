// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <span>
#include <vector>

namespace Acts {

class Surface;

/// The NavigationStream is a container for the navigation candidates that
/// are currentlu processed in a given context. The context could be local to a
/// volume, or global to an entire track following.
///
/// The current candidates are stored in a vector of candidates, where an index
/// is used to indicate the current active candidate.
class NavigationStream {
 public:
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

  /// Switch to next next candidate
  ///
  /// @return true if a next candidate is available
  bool switchToNextCandidate() {
    if (m_currentIndex < m_candidates.size()) {
      ++m_currentIndex;
      return true;
    }
    return false;
  }

  /// Const access the current candidate
  /// @return Const reference to current candidate
  const NavigationTarget& currentCandidate() const {
    return m_candidates.at(m_currentIndex);
  }

  /// Current Index
  /// @return Index of the current candidate in the vector
  std::size_t currentIndex() const { return m_currentIndex; }

  /// Non-const access the candidate vector
  /// @return Mutable reference to vector of navigation candidates
  std::vector<NavigationTarget>& candidates() { return m_candidates; }

  /// Const access the candidate vector
  /// @return Const reference to vector of navigation candidates
  const std::vector<NavigationTarget>& candidates() const {
    return m_candidates;
  }

  /// Non-const access the current candidate
  ///
  /// This will throw and out of bounds exception if the stream is not
  /// valid anymore.
  /// @return Mutable reference to current candidate
  NavigationTarget& currentCandidate() {
    return m_candidates.at(m_currentIndex);
  }

  /// The number of active candidates
  /// @return Number of remaining candidates from current position onwards
  std::size_t remainingCandidates() const {
    return (m_candidates.size() - m_currentIndex);
  }

  /// Fill one surface into the candidate vector
  ///
  /// @param surface the surface to be filled
  /// @param bTolerance the boundary tolerance used for the intersection
  void addSurfaceCandidate(const Surface& surface,
                           const BoundaryTolerance& bTolerance);

  /// Fill n surfaces into the candidate vector
  ///
  /// @param surfaces the surfaces that are filled in
  /// @param bTolerance the boundary tolerance used for the intersection
  void addSurfaceCandidates(std::span<const Surface*> surfaces,
                            const BoundaryTolerance& bTolerance);

  /// Fill one portal into the candidate vector
  ///
  /// @param portal the portals that are filled in
  void addPortalCandidate(const Portal& portal);

  /// Initialize the stream from a query point
  ///
  /// @param gctx is the geometry context
  /// @param queryPoint holds current position, direction, etc.
  /// @param cTolerance is the candidate search tolerance
  /// @param onSurfaceTolerance is the tolerance for on-surface intersections
  ///
  /// This method will first de-duplicate the candidates on basis of the surface
  /// pointer to make sure that the multi-intersections are handled correctly.
  /// This will allow intializeStream() to be called even as a re-initialization
  /// and still work correctly with at one time valid candidates.
  ///
  /// @return true if the stream is active, false indicates that there are no valid candidates
  bool initialize(const GeometryContext& gctx,
                  const NavigationStream::QueryPoint& queryPoint,
                  const BoundaryTolerance& cTolerance,
                  double onSurfaceTolerance = s_onSurfaceTolerance);

  /// Convenience method to update a stream from a new query point,
  /// this could be called from navigation delegates that do not require
  /// a local state or from the navigator on the target stream
  ///
  /// @param gctx is the geometry context
  /// @param queryPoint holds current position, direction, etc.
  /// @param onSurfaceTolerance is the tolerance for on-surface intersections
  ///
  /// @return true if the stream is active, false indicate no valid candidates left
  bool update(const GeometryContext& gctx,
              const NavigationStream::QueryPoint& queryPoint,
              double onSurfaceTolerance = s_onSurfaceTolerance);

  /// Reset the navigation stream by clearing all candidates and resetting the
  /// index.
  ///
  /// This clears the candidates vector and resets the current index to 0.
  void reset();

 private:
  /// The candidates of this navigation stream
  std::vector<NavigationTarget> m_candidates;

  /// The currently active candidate
  std::size_t m_currentIndex = 0u;
};

/// Append-only helper to add candidates to a navigation stream.
class AppendOnlyNavigationStream {
 public:
  /// Constructor from navigation stream reference
  /// @param stream Navigation stream to append to
  explicit AppendOnlyNavigationStream(NavigationStream& stream);
  /// Add a surface candidate to the stream
  /// @param surface The surface to add
  /// @param bTolerance Boundary tolerance for the surface
  void addSurfaceCandidate(const Surface& surface,
                           const BoundaryTolerance& bTolerance);
  /// Add a portal candidate to the stream
  /// @param portal The portal to add
  void addPortalCandidate(const Portal& portal);

 private:
  NavigationStream* m_stream;
};

}  // namespace Acts
