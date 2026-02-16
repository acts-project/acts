// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include <map>
#include <optional>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief The navigation options for the tracking geometry
///
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// Hint for start object
  const object_t* startObject = nullptr;
  /// Hint for end object
  const object_t* endObject = nullptr;

  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The minimum distance for a surface to be considered
  double nearLimit = 0;
  /// The maximum distance for a surface to be considered
  double farLimit = std::numeric_limits<double>::max();
};

/// @brief Steers the propagation through the geometry by providing the next
/// surface to be targeted.
///
/// The Navigator is part of the propagation and responsible for steering
/// the surface sequence to encounter all the relevant surfaces which are
/// intersected by the trajectory.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation flow is hit, it is
/// set to the navigation state, such that other actors can deal with it.
///
/// The current target surface is referenced by an index which points into
/// the navigation candidates. The navigation candidates are ordered by the
/// path length to the surface. If a surface is hit, the
/// `state.currentSurface` pointer is set. This actors to observe
/// that we are on a surface.
///
class Navigator {
 public:
  /// Type alias for navigation surface candidates container
  using NavigationSurfaces =
      boost::container::small_vector<NavigationTarget, 10>;

  /// Type alias for navigation layer candidates container
  using NavigationLayers = boost::container::small_vector<NavigationTarget, 10>;

  /// Type alias for navigation boundary candidates container
  using NavigationBoundaries =
      boost::container::small_vector<NavigationTarget, 4>;

  /// Type alias for generic navigation candidates container
  using NavigationCandidates =
      boost::container::small_vector<NavigationTarget, 10>;

  /// Type alias for external surfaces container
  using ExternalSurfaces = std::vector<GeometryIdentifier>;

  /// Type alias for geometry version enumeration
  using GeometryVersion = TrackingGeometry::GeometryVersion;

  /// The navigation stage
  enum struct Stage : int {
    initial = 0,
    surfaceTarget = 1,
    layerTarget = 2,
    boundaryTarget = 3,
  };

  /// The navigator configuration
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;
  };

  /// The navigator options
  struct Options : public NavigatorPlainOptions {
    /// Constructor with geometry context
    /// @param gctx The geometry context for the navigation
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};
    /// Surfaces that are not part of the tracking geometry
    std::vector<const Surface*> freeSurfaces = {};

    /// Insert an external surface to be considered during navigation
    /// @param surface: The surface to add to the list
    void insertExternalSurface(const Surface& surface) {
      if (surface.geometryId() != GeometryIdentifier{}) {
        externalSurfaces.push_back(surface.geometryId());
      } else {
        freeSurfaces.push_back(&surface);
      }
    }
    /// @brief Delegate to decide whether free surfaces are appended to the navigation
    ///        stream given the current volume and the track coordinates. If the
    ///        delegate is set, it is called in each candidate resolution step
    ///        for each surface that has not been marked as reached yet.
    /// @param gctx: Current geometry context carrying the alignment information
    /// @param currentVol: The current tracking volume in which the propagator resides
    /// @param pos: Position of the track in global coordinates
    /// @param dir: Direction vector of the track
    /// @param surface: Free surface candidate to test
    using FreeSurfaceSelctor_t = Delegate<bool(
        const GeometryContext& gctx, const TrackingVolume& currentVol,
        const Vector3& pos, const Vector3& dir, const Surface& candidate)>;
    /// Delegate for selecting free surfaces during navigation
    FreeSurfaceSelctor_t freeSurfaceSelector{};
    /// Set the plain navigation options
    /// @param options The plain navigator options to set
    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    /// Constructor with navigation options
    /// @param options_ The navigation options for this state
    explicit State(const Options& options_) : options(options_) {}

    /// Navigation options configuration
    Options options;

    /// Management of policy state allocation and deallocation
    NavigationPolicyStateManager policyStateManager;

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface index of the navigation state
    std::optional<std::size_t> navSurfaceIndex;

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer index of the navigation state
    std::optional<std::size_t> navLayerIndex;

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary index of the navigation state
    std::optional<std::size_t> navBoundaryIndex;

    // Navigation candidates(portals and surfaces together)
    /// the vector of navigation candidates to work through
    NavigationCandidates navCandidates = {};
    /// the current candidate index of the navigation state
    std::optional<std::size_t> navCandidateIndex;

    /// Free candidates not part of the tracking geometry.
    //  They are stored as a pair of surface pointer
    /// and a boolean indicating whether the surface has already been
    /// reached during propagation
    std::vector<std::pair<const Surface*, bool>> freeCandidates{};

    /// Get reference to current navigation surface
    /// @return Reference to current navigation target
    NavigationTarget& navSurface() {
      return navSurfaces.at(navSurfaceIndex.value());
    }

    /// Get reference to current navigation layer
    /// @return Reference to current layer intersection
    NavigationTarget& navLayer() { return navLayers.at(navLayerIndex.value()); }

    /// Get reference to current navigation boundary
    /// @return Reference to current boundary intersection
    NavigationTarget& navBoundary() {
      return navBoundaries.at(navBoundaryIndex.value());
    }

    /// Get reference to current navigation candidate
    /// @return Reference to current boundary intersection
    NavigationTarget& navCandidate() {
      return navCandidates.at(navCandidateIndex.value());
    }

    /// Volume where the navigation started
    const TrackingVolume* startVolume = nullptr;
    /// Layer where the navigation started
    const Layer* startLayer = nullptr;
    /// Surface where the navigation started
    const Surface* startSurface = nullptr;
    /// Current volume during navigation
    const TrackingVolume* currentVolume = nullptr;
    /// Current layer during navigation
    const Layer* currentLayer = nullptr;
    /// Current surface during navigation
    const Surface* currentSurface = nullptr;
    /// Target surface for navigation
    const Surface* targetSurface = nullptr;

    /// Flag to break navigation loop
    bool navigationBreak = false;
    /// Current navigation stage in the state machine
    Stage navigationStage = Stage::initial;

    /// Statistics collection for navigation performance
    NavigatorStatistics statistics;

    /// Stream for navigation debugging and monitoring
    NavigationStream stream;

    /// Reset navigation state after switching layers
    void resetAfterLayerSwitch() {
      navSurfaces.clear();
      navSurfaceIndex.reset();
    }

    /// Reset navigation state after switching volumes
    void resetAfterVolumeSwitch() {
      resetAfterLayerSwitch();

      navLayers.clear();
      navLayerIndex.reset();
      navBoundaries.clear();
      navBoundaryIndex.reset();
      navCandidates.clear();
      navCandidateIndex.reset();

      currentLayer = nullptr;

      policyStateManager.reset();
    }

    /// Completely reset navigation state to initial conditions
    void reset() {
      resetAfterVolumeSwitch();

      currentVolume = nullptr;
      currentSurface = nullptr;

      navigationBreak = false;
      navigationStage = Stage::initial;
      // Set the surface reached switches back to false
      std::ranges::for_each(freeCandidates,
                            [](std::pair<const Surface*, bool>& freeSurface) {
                              freeSurface.second = false;
                            });
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO));

  /// Create a navigation state from options
  /// @param options The navigation options
  /// @return A new navigation state
  State makeState(const Options& options) const;

  /// Get the current surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to current surface, or nullptr if none
  const Surface* currentSurface(const State& state) const;

  /// Get the current volume from navigation state
  /// @param state The navigation state
  /// @return Pointer to current volume, or nullptr if none
  const TrackingVolume* currentVolume(const State& state) const;

  /// Get material properties of the current volume
  /// @param state The navigation state
  /// @return Pointer to volume material, or nullptr if no volume or material
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const;

  /// Get the starting surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to start surface, or nullptr if none
  const Surface* startSurface(const State& state) const;

  /// Get the target surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to target surface, or nullptr if none
  const Surface* targetSurface(const State& state) const;

  /// Check if navigation has reached the end of the world (no current volume)
  /// @param state The navigation state
  /// @return True if end of world is reached
  bool endOfWorldReached(const State& state) const;

  /// Check if navigation should be interrupted
  /// @param state The navigation state
  /// @return True if navigation break flag is set
  bool navigationBreak(const State& state) const;

  /// @brief Initialize the navigator state
  ///
  /// This function initializes the navigator state for a new propagation.
  ///
  /// @param state The navigation state
  /// @param position The start position
  /// @param direction The start direction
  /// @param propagationDirection The propagation direction
  ///
  /// @return Indication if the initialization was successful
  [[nodiscard]] Result<void> initialize(State& state, const Vector3& position,
                                        const Vector3& direction,
                                        Direction propagationDirection) const;

  /// @brief Get the next target surface
  ///
  /// This function gets the next target surface for the propagation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const;

  /// @brief Check if the current target is still valid
  ///
  /// This function checks if the target is valid.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is valid
  bool checkTargetValid(State& state, const Vector3& position,
                        const Vector3& direction) const;

  /// @brief Handle the surface reached
  ///
  /// This function handles the surface reached.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param surface The surface reached
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const;

 private:
  /// @brief NextTarget helper function for Gen1 geometry configuration
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  NavigationTarget getNextTargetGen1(State& state, const Vector3& position,
                                     const Vector3& direction) const;

  /// @brief NextTarget helper function for Gen3 geometry configuration
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  NavigationTarget getNextTargetGen3(State& state, const Vector3& position,
                                     const Vector3& direction) const;

  /// @brief NextTarget helper function
  /// This function is called for returning the next target
  /// and checks gen1/gen3 case in order to sub-call the proper functions
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  NavigationTarget tryGetNextTarget(State& state, const Vector3& position,
                                    const Vector3& direction) const;

  /// @brief Resolve compatible candidates (surfaces or portals) for gen3
  /// navigation
  ///
  /// This function is called when gen3 configuration is found and it resolves
  /// at the same time for portals and surfaces
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveCandidates(State& state, const Vector3& position,
                         const Vector3& direction) const;

  /// @brief Resolve compatible surfaces
  ///
  /// This function resolves the compatible surfaces for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveSurfaces(State& state, const Vector3& position,
                       const Vector3& direction) const;

  /// @brief Resolve compatible layers
  ///
  /// This function resolves the compatible layers for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveLayers(State& state, const Vector3& position,
                     const Vector3& direction) const;

  /// @brief Resolve compatible boundaries
  ///
  /// This function resolves the compatible boundaries for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveBoundaries(State& state, const Vector3& position,
                         const Vector3& direction) const;

  /// @brief Check if the navigator is inactive
  ///
  /// This function checks if the navigator is inactive.
  ///
  /// @param state The navigation state
  ///
  /// @return True if the navigator is inactive
  bool inactive(const State& state) const;

  /// @brief Get volume info string for logging
  ///
  /// @tparam propagator_state_t The propagator state type
  /// @param state The state containing current volume info
  /// @return String with volume name for logging
  std::string volInfo(const State& state) const;

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  // Cached so we don't have to query the TrackingGeometry constantly.
  TrackingGeometry::GeometryVersion m_geometryVersion{};

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
