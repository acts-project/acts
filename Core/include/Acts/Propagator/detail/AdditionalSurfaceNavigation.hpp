// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <array>
#include <stdexcept>
#include <vector>

namespace Acts::detail {

/// Offers the additional surfaces of the navigator options and the target
/// surface, on top of the candidates of a geometry navigator. Each navigation
/// state owns one instance.
///
/// The additional surfaces and the candidates of the geometry navigator are two
/// streams. Both are intersected again before every step, and the closer one is
/// the target. An intersection seen ahead stays a candidate when it moves
/// behind, so the propagation can correct backwards until the surface is
/// reached. Unreachable targets are bounded by the propagation limits.
class AdditionalSurfaceNavigation {
  /// An additional surface with its bookkeeping
  struct AdditionalSurfaceState {
    /// The additional surface
    AdditionalSurface surface;
    /// Whether the propagation reached the surface
    bool reached = false;
    /// Once seen ahead, keep each intersection available until the surface
    /// is reached, even if it moves behind the propagation.
    std::array<bool, s_maximumNumberOfIntersections> pending{};
  };

  /// Additional surfaces of the options and the target surface
  std::vector<AdditionalSurfaceState> m_additionalSurfaces;

  /// The candidate of another navigator an additional surface took
  /// precedence over. It is handed out once the additional surface is no
  /// longer the closer one.
  NavigationTarget m_heldBack = NavigationTarget::None();

  /// The target handed out to the propagation
  NavigationTarget m_handedOut = NavigationTarget::None();
  /// Whether the target handed out is an additional surface
  bool m_handedOutIsAdditional = false;

 public:
  /// Reset the surface list and candidate history for a navigation run.
  /// @param options The navigator options
  /// @param targetSurface Optional target to include with its bounds checked
  void initialize(const NavigatorPlainOptions& options,
                  const Surface* targetSurface) {
    m_additionalSurfaces.clear();
    m_additionalSurfaces.reserve(options.additionalSurfaces.size() + 1);
    for (const AdditionalSurface& surface : options.additionalSurfaces) {
      if (surface.surface == nullptr) {
        throw std::invalid_argument("Navigator: additional surface is nullptr");
      }
      m_additionalSurfaces.push_back({surface});
    }
    // Propagation targets require a bounds check.
    if (targetSurface != nullptr) {
      m_additionalSurfaces.push_back(
          {{targetSurface, BoundaryTolerance::None()}});
    }

    m_heldBack = NavigationTarget::None();
    m_handedOut = NavigationTarget::None();
    m_handedOutIsAdditional = false;
  }

  /// Returns the closest additional surface
  /// @param options The navigator options
  /// @param position The current position
  /// @param direction The current direction
  /// @return The next target, or none
  NavigationTarget nextTarget(const NavigatorPlainOptions& options,
                              const Vector3& position,
                              const Vector3& direction) {
    return handOut(closest(options, position, direction, nullptr), true);
  }

  /// Returns the closer of the next candidate of another navigator and the
  /// closest additional surface. The candidate is held back, not skipped,
  /// while an additional surface is closer.
  ///
  /// @param options The navigator options
  /// @param position The current position
  /// @param direction The current direction
  /// @param volume The current volume of the other navigator, read after
  ///        @p nextCandidate ran
  /// @param nextCandidate Callable returning the next candidate of the other
  ///        navigator
  /// @return The next target
  template <typename next_candidate_t>
  NavigationTarget nextTarget(const NavigatorPlainOptions& options,
                              const Vector3& position, const Vector3& direction,
                              const TrackingVolume* const& volume,
                              next_candidate_t&& nextCandidate) {
    if (m_additionalSurfaces.empty()) {
      return nextCandidate();
    }

    // Ask the other navigator again if it had no candidate
    if (m_heldBack.isNone()) {
      m_heldBack = nextCandidate();
    }

    if (const NavigationTarget additional =
            closest(options, position, direction, volume);
        !additional.isNone()) {
      refresh(options, m_heldBack, position, direction);
      // A tie goes to the candidate of the other navigator
      if (m_heldBack.isNone() ||
          additional.pathLength() < m_heldBack.intersection().pathLength()) {
        return handOut(additional, true);
      }
    }

    const NavigationTarget candidate = m_heldBack;
    m_heldBack = NavigationTarget::None();
    return handOut(candidate, false);
  }

  /// Rechecks the selected target against both candidate streams. Pending
  /// additional surfaces remain candidates when their intersections move
  /// behind.
  ///
  /// @param options The navigator options
  /// @param position The current position
  /// @param direction The current direction
  /// @param volume The current volume of the other navigator
  /// @param candidateValid Whether the other navigator holds its candidates
  ///        valid. If not, it is asked for a new one.
  /// @return True if the target handed out is still the closest
  bool checkTargetValid(const NavigatorPlainOptions& options,
                        const Vector3& position, const Vector3& direction,
                        const TrackingVolume* volume, bool candidateValid) {
    if (!candidateValid) {
      m_heldBack = NavigationTarget::None();
      handOut(NavigationTarget::None(), false);
      return false;
    }
    if (m_additionalSurfaces.empty() || m_handedOut.isNone()) {
      return true;
    }

    const NavigationTarget additional =
        closest(options, position, direction, volume);

    if (m_handedOutIsAdditional) {
      if (additional.isNone() ||
          &additional.surface() != &m_handedOut.surface() ||
          additional.intersectionIndex() != m_handedOut.intersectionIndex()) {
        return false;
      }
      if (!m_heldBack.isNone()) {
        refresh(options, m_heldBack, position, direction);
        // A tie goes to the candidate of the other navigator
        if (m_heldBack.intersection().pathLength() <= additional.pathLength()) {
          return false;
        }
      }
      return true;
    }

    if (additional.isNone()) {
      return true;
    }
    refresh(options, m_handedOut, position, direction);
    if (additional.pathLength() < m_handedOut.intersection().pathLength()) {
      // Hold the candidate back instead of asking the other navigator again,
      // which would move on to its next candidate
      m_heldBack = m_handedOut;
      m_handedOut = NavigationTarget::None();
      return false;
    }
    return true;
  }

  /// Handles reaching a surface for another navigator
  /// @param surface The surface reached
  /// @return True if the surface is an additional surface the other navigator
  ///         did not target, so the other navigator has nothing to do
  bool handleSurfaceReached(const Surface& surface) {
    bool found = false;
    for (AdditionalSurfaceState& additional : m_additionalSurfaces) {
      if (additional.surface.surface == &surface) {
        additional.reached = true;
        additional.pending.fill(false);
        found = true;
      }
    }
    return found && m_handedOutIsAdditional;
  }

 private:
  /// Records the target handed out
  NavigationTarget handOut(const NavigationTarget& target, bool isAdditional) {
    m_handedOut = target;
    m_handedOutIsAdditional = isAdditional && !target.isNone();
    return target;
  }

  /// Intersects a candidate of the other navigator again. Keeps the stored
  /// intersection if the straight-line estimate lost the surface.
  static void refresh(const NavigatorPlainOptions& options,
                      NavigationTarget& candidate, const Vector3& position,
                      const Vector3& direction) {
    if (candidate.isNone()) {
      return;
    }
    const Intersection3D refreshed =
        candidate.surface()
            .intersect(options.geoContext, position, direction,
                       candidate.boundaryTolerance(), options.surfaceTolerance)
            .at(candidate.intersectionIndex());
    if (refreshed.isValid()) {
      candidate.intersection() = refreshed;
    }
  }

  /// Closest additional surface in the current volume, ahead or pending.
  /// Intersections are tracked by index, which for a surface with two
  /// solutions is their order in signed path length.
  NavigationTarget closest(const NavigatorPlainOptions& options,
                           const Vector3& position, const Vector3& direction,
                           const TrackingVolume* volume) {
    NavigationTarget closest = NavigationTarget::None();

    for (AdditionalSurfaceState& additional : m_additionalSurfaces) {
      const AdditionalSurface& surface = additional.surface;
      if (additional.reached && surface.dropAfterReached) {
        continue;
      }
      if (surface.volume != nullptr && surface.volume != volume) {
        continue;
      }

      // The closer solution can be behind the propagation, so check all
      const MultiIntersection3D multiIntersection = surface.surface->intersect(
          options.geoContext, position, direction, surface.boundaryTolerance,
          options.surfaceTolerance);

      for (const auto [intersectionIndex, intersection] :
           enumerate(multiIntersection)) {
        if (!intersection.isValid() ||
            intersection.pathLength() > options.farLimit) {
          continue;
        }
        const double pathLength = intersection.pathLength();
        const bool onSurface =
            intersection.status() == IntersectionStatus::onSurface;
        if (!onSurface && pathLength > 0) {
          additional.pending[intersectionIndex] = true;
        }

        // Retain an earlier forward candidate even if it is now behind, as a
        // turn can change the order of the candidates. The stepper corrects
        // backwards, just as for a direct navigation target. A kept surface
        // is not reached twice at the same position.
        const bool offer = onSurface ? !additional.reached
                                     : additional.pending[intersectionIndex];
        if (!offer) {
          continue;
        }
        if (closest.isNone() || pathLength < closest.pathLength()) {
          closest = NavigationTarget(
              intersection, static_cast<IntersectionIndex>(intersectionIndex),
              *surface.surface, surface.boundaryTolerance);
        }
      }
    }

    return closest;
  }
};

}  // namespace Acts::detail
