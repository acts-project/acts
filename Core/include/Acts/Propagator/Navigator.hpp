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
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorError.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <map>
#include <optional>
#include <sstream>
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

  /// Type alias for external surfaces map indexed by layer ID
  using ExternalSurfaces = std::multimap<std::uint64_t, GeometryIdentifier>;

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

    /// Insert an external surface to be considered during navigation
    /// @param geoid Geometry identifier of the surface to insert
    void insertExternalSurface(GeometryIdentifier geoid) {
      externalSurfaces.insert(
          std::pair<std::uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
    }

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

    /// Get reference to current navigation surface
    /// @return Reference to current surface target
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

    /// Gen3 surface stream
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

      currentLayer = nullptr;
    }

    /// Completely reset navigation state to initial conditions
    void reset() {
      resetAfterVolumeSwitch();

      currentVolume = nullptr;
      currentSurface = nullptr;

      navigationBreak = false;
      navigationStage = Stage::initial;

      stream.reset();
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {
    if (m_cfg.trackingGeometry == nullptr) {
      throw std::invalid_argument("Navigator: No tracking geometry provided.");
    }
    m_geometryVersion = m_cfg.trackingGeometry->geometryVersion();
  }

  /// Create a navigation state from options
  /// @param options The navigation options
  /// @return A new navigation state
  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  /// Get the current surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to current surface, or nullptr if none
  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  /// Get the current volume from navigation state
  /// @param state The navigation state
  /// @return Pointer to current volume, or nullptr if none
  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  /// Get material properties of the current volume
  /// @param state The navigation state
  /// @return Pointer to volume material, or nullptr if no volume or material
  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  /// Get the starting surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to start surface, or nullptr if none
  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  /// Get the target surface from navigation state
  /// @param state The navigation state
  /// @return Pointer to target surface, or nullptr if none
  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  /// Check if navigation has reached the end of the world (no current volume)
  /// @param state The navigation state
  /// @return True if end of world is reached
  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  /// Check if navigation should be interrupted
  /// @param state The navigation state
  /// @return True if navigation break flag is set
  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

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
                                        Direction propagationDirection) const {
    (void)propagationDirection;

    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    auto printGeometryVersion = [](auto ver) {
      using enum TrackingGeometry::GeometryVersion;
      switch (ver) {
        case Gen1:
          return "Gen1";
        case Gen3:
          return "Gen3";
        default:
          throw std::runtime_error("Unknown geometry version.");
      }
    };
    ACTS_VERBOSE(volInfo(state) << "Geometry version is: "
                                << printGeometryVersion(m_geometryVersion));

    state.reset();

    if (m_geometryVersion == GeometryVersion::Gen3) {
      // Empirical pre-allocation of candidates for the next navigation
      // iteration.
      // @TODO: Make this user configurable through the configuration
      state.stream.candidates().reserve(50);
    }

    state.startSurface = state.options.startSurface;
    state.targetSurface = state.options.targetSurface;

    // @TODO: Implement fast initialization with Gen3. This requires the volume lookup to work properly

    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Surface.");

      state.startLayer = state.startSurface->associatedLayer();
      state.startVolume = state.startLayer->trackingVolume();
    } else if (state.startVolume != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Volume.");

      state.startLayer = state.startVolume->associatedLayer(
          state.options.geoContext, position);
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Slow start initialization through search.");
      ACTS_VERBOSE(volInfo(state)
                   << "Starting from position " << toString(position)
                   << " and direction " << toString(direction));

      // current volume and layer search through global search
      state.startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);

      if (state.startVolume != nullptr) {
        state.startLayer = state.startVolume->associatedLayer(
            state.options.geoContext, position);
      } else {
        ACTS_ERROR(volInfo(state)
                   << "No start volume resolved. Nothing left to do.");
        state.navigationBreak = true;
      }
    }

    state.currentVolume = state.startVolume;
    state.currentLayer = state.startLayer;
    state.currentSurface = state.startSurface;

    if (state.currentVolume != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start volume resolved "
                                  << state.currentVolume->geometryId());

      if (!state.currentVolume->inside(position,
                                       state.options.surfaceTolerance)) {
        ACTS_DEBUG(
            volInfo(state)
            << "We did not end up inside the expected volume. position = "
            << position.transpose());

        return Result<void>::failure(NavigatorError::NotInsideExpectedVolume);
      }
    }
    if (state.currentLayer != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start layer resolved "
                                  << state.currentLayer->geometryId());
    }
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start surface resolved "
                                  << state.currentSurface->geometryId());

      if (!state.currentSurface->isOnSurface(
              state.options.geoContext, position, direction,
              BoundaryTolerance::Infinite(), state.options.surfaceTolerance)) {
        ACTS_DEBUG(volInfo(state)
                   << "We did not end up on the expected surface. surface = "
                   << state.currentSurface->geometryId()
                   << " position = " << position.transpose()
                   << " direction = " << direction.transpose());

        return Result<void>::failure(NavigatorError::NotOnExpectedSurface);
      }
    }

    return Result<void>::success();
  }

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
                              const Vector3& direction) const {
    // Reset the current surface
    state.currentSurface = nullptr;

    if (inactive(state)) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::nextTarget.");

    auto tryGetNextTarget = [&]() -> NavigationTarget {
      // Try targeting the surfaces - then layers - then boundaries

      if (state.navigationStage == Stage::initial) {
        ACTS_VERBOSE(volInfo(state) << "Target surfaces.");
        state.navigationStage = Stage::surfaceTarget;
      }

      if (state.navigationStage == Stage::surfaceTarget) {
        if (!state.navSurfaceIndex.has_value()) {
          // First time, resolve the surfaces
          resolveSurfaces(state, position, direction);
          state.navSurfaceIndex = 0;
        } else {
          ++state.navSurfaceIndex.value();
        }
        if (state.navSurfaceIndex.value() < state.navSurfaces.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
          return state.navSurface();
        } else {
          // This was the last surface, switch to layers
          ACTS_VERBOSE(volInfo(state) << "Target layers.");
          if (m_geometryVersion == GeometryVersion::Gen1) {
            state.navigationStage = Stage::layerTarget;
          } else {
            state.navigationStage = Stage::boundaryTarget;
          }
        }
      }

      if (state.navigationStage == Stage::layerTarget) {
        if (!state.navLayerIndex.has_value()) {
          // First time, resolve the layers
          resolveLayers(state, position, direction);
          state.navLayerIndex = 0;
        } else {
          ++state.navLayerIndex.value();
        }
        if (state.navLayerIndex.value() < state.navLayers.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
          return state.navLayer();
        } else {
          // This was the last layer, switch to boundaries
          ACTS_VERBOSE(volInfo(state) << "Target boundaries.");
          state.navigationStage = Stage::boundaryTarget;
        }
      }

      if (state.navigationStage == Stage::boundaryTarget) {
        if (!state.navBoundaryIndex.has_value()) {
          // First time, resolve the boundaries
          resolveBoundaries(state, position, direction);
          state.navBoundaryIndex = 0;
        } else {
          ++state.navBoundaryIndex.value();
        }
        if (state.navBoundaryIndex.value() < state.navBoundaries.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
          return state.navBoundary();
        } else {
          // This was the last boundary, we have to leave the volume somehow,
          // renavigate
          ACTS_VERBOSE(volInfo(state)
                       << "Boundary targets exhausted. Renavigate.");
        }
      }

      ACTS_VERBOSE(volInfo(state)
                   << "Unknown state. No target found. Renavigate.");
      return NavigationTarget::None();
    };

    NavigationTarget nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    state.reset();
    ++state.statistics.nRenavigations;

    // We might have punched through a boundary and entered another volume
    // so we have to reinitialize
    state.currentVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
        state.options.geoContext, position);

    if (state.currentVolume == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "No volume found, stop navigation.");
      state.navigationBreak = true;
      return NavigationTarget::None();
    }

    state.currentLayer = state.currentVolume->associatedLayer(
        state.options.geoContext, position);

    ACTS_VERBOSE(volInfo(state) << "Resolved volume and layer.");

    // Rerun the targeting
    nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    ACTS_VERBOSE(
        volInfo(state)
        << "No targets found again, we got really lost! Stop navigation.");
    state.navigationBreak = true;
    return NavigationTarget::None();
  }

  /// @brief Check if the current target is still valid
  ///
  /// This function checks if the target is valid.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)position;
    (void)direction;

    return state.navigationStage != Stage::initial;
  }

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
                            const Surface& surface) const {
    if (inactive(state)) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::handleSurfaceReached.");

    state.currentSurface = &surface;

    ACTS_VERBOSE(volInfo(state)
                 << "Current surface: " << state.currentSurface->geometryId());

    if (state.navigationStage == Stage::surfaceTarget &&
        &state.navSurface().surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling surface status.");

      return;
    }

    if (state.navigationStage == Stage::layerTarget &&
        &state.navLayer().surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling layer status.");

      // Switch to the next layer
      state.currentLayer = &state.navLayer().layer();
      state.navigationStage = Stage::surfaceTarget;

      // partial reset
      state.resetAfterLayerSwitch();

      return;
    }

    if (state.navigationStage == Stage::boundaryTarget &&
        &state.navBoundary().surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling boundary status.");

      if (m_geometryVersion == GeometryVersion::Gen1) {
        // Switch to the next volume using the boundary
        const BoundarySurface* boundary =
            &state.navBoundary().boundarySurface();
        assert(boundary != nullptr && "Retrieved boundary surface is nullptr");
        state.currentVolume = boundary->attachedVolume(state.options.geoContext,
                                                       position, direction);
      } else {
        const Portal* portal = &state.navBoundary().portal();
        assert(portal != nullptr && "Retrieved portal is nullptr");
        auto res = portal->resolveVolume(state.options.geoContext, position,
                                         direction);
        if (!res.ok()) {
          ACTS_ERROR(volInfo(state)
                     << "Failed to resolve volume through portal: "
                     << res.error().message());
          return;
        }

        state.currentVolume = res.value();
      }

      // partial reset
      state.resetAfterVolumeSwitch();

      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Volume updated.");
        if (m_geometryVersion == GeometryVersion::Gen1) {
          state.navigationStage = Stage::layerTarget;
        } else {
          state.navigationStage = Stage::surfaceTarget;
        }
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "No more volume to progress to, stopping navigation.");
        state.navigationBreak = true;
      }

      return;
    }

    ACTS_ERROR(volInfo(state) << "Surface reached but unknown state.");
  }

 private:
  /// @brief Resolve compatible surfaces
  ///
  /// This function resolves the compatible surfaces for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveSurfaces(State& state, const Vector3& position,
                       const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible surfaces.");

    if (m_geometryVersion == GeometryVersion::Gen1) {
      const Layer* currentLayer = state.currentLayer;

      if (currentLayer == nullptr) {
        ACTS_VERBOSE(volInfo(state) << "No layer to resolve surfaces.");
        return;
      }

      const Surface* layerSurface = &currentLayer->surfaceRepresentation();

      NavigationOptions<Surface> navOpts;
      navOpts.resolveSensitive = m_cfg.resolveSensitive;
      navOpts.resolveMaterial = m_cfg.resolveMaterial;
      navOpts.resolvePassive = m_cfg.resolvePassive;
      navOpts.startObject = state.currentSurface;
      navOpts.endObject = state.targetSurface;
      navOpts.nearLimit = state.options.nearLimit;
      navOpts.farLimit = state.options.farLimit;

      if (!state.options.externalSurfaces.empty()) {
        auto layerId = layerSurface->geometryId().layer();
        auto externalSurfaceRange =
            state.options.externalSurfaces.equal_range(layerId);
        navOpts.externalSurfaces.reserve(
            state.options.externalSurfaces.count(layerId));
        for (auto itSurface = externalSurfaceRange.first;
             itSurface != externalSurfaceRange.second; itSurface++) {
          navOpts.externalSurfaces.push_back(itSurface->second);
        }
      }

      // Request the compatible surfaces
      state.navSurfaces = currentLayer->compatibleSurfaces(
          state.options.geoContext, position, direction, navOpts);
      // Sort the surfaces by path length.
      // Special care is taken for the external surfaces which should always
      // come first, so they are preferred to be targeted and hit first.
      std::ranges::sort(state.navSurfaces, [&state](const NavigationTarget& a,
                                                    const NavigationTarget& b) {
        // Prefer to sort by path length. We assume surfaces are at the same
        // distance if the difference is smaller than the tolerance.
        if (std::abs(a.pathLength() - b.pathLength()) >
            state.options.surfaceTolerance) {
          return NavigationTarget::pathLengthOrder(a, b);
        }
        // If the path length is practically the same, sort by geometry.
        // First we check if one of the surfaces is external.
        bool aIsExternal = a.boundaryTolerance().isInfinite();
        bool bIsExternal = b.boundaryTolerance().isInfinite();
        if (aIsExternal == bIsExternal) {
          // If both are external or both are not external, sort by geometry
          // identifier
          return a.surface().geometryId() < b.surface().geometryId();
        }
        // If only one is external, it should come first
        return aIsExternal;
      });
      // For now we implicitly remove overlapping surfaces.
      // For track finding it might be useful to discover overlapping surfaces
      // and check for compatible measurements. This is under investigation
      // and might be implemented in the future.
      auto toBeRemoved = std::ranges::unique(
          state.navSurfaces, [&](const auto& a, const auto& b) {
            return std::abs(a.pathLength() - b.pathLength()) <
                   state.options.surfaceTolerance;
          });
      if (toBeRemoved.begin() != toBeRemoved.end()) {
        ACTS_VERBOSE(volInfo(state)
                     << "Removing "
                     << std::distance(toBeRemoved.begin(), toBeRemoved.end())
                     << " overlapping surfaces.");
      }
      state.navSurfaces.erase(toBeRemoved.begin(), toBeRemoved.end());
    } else {
      // @TODO: What to do with external surfaces?
      // Gen 3 !
      state.stream.reset();
      AppendOnlyNavigationStream appendOnly{state.stream};
      NavigationArguments args;
      args.position = position;
      args.direction = direction;
      args.wantsPortals = false;
      args.wantsSurfaces = true;
      state.currentVolume->initializeNavigationCandidates(args, appendOnly,
                                                          logger());

      // Filter out portals before intersection

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << state.stream.candidates().size()
                   << " navigation candidates.");

      state.stream.initialize(state.options.geoContext, {position, direction},
                              BoundaryTolerance::None(),
                              state.options.surfaceTolerance);
      ACTS_VERBOSE(volInfo(state)
                   << "Now " << state.stream.candidates().size()
                   << " navigation candidates after initialization");

      state.navSurfaces.clear();

      auto it = std::ranges::find_if(
          state.stream.candidates(), [&](const auto& candidate) {
            return detail::checkPathLength(candidate.pathLength(),
                                           state.options.nearLimit,
                                           state.options.farLimit, logger());
          });

      std::ranges::copy(
          std::ranges::subrange(it, state.stream.candidates().end()),
          std::back_inserter(state.navSurfaces));
    }

    // Print surface information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navSurfaces.size();
      os << " surface candidates found at path(s): ";
      for (auto& sfc : state.navSurfaces) {
        os << sfc.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navSurfaces.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
    }
  }

  /// @brief Resolve compatible layers
  ///
  /// This function resolves the compatible layers for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveLayers(State& state, const Vector3& position,
                     const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

    NavigationOptions<Layer> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.currentLayer;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    // Request the compatible layers
    state.navLayers = state.currentVolume->compatibleLayers(
        state.options.geoContext, position, direction, navOpts);
    std::ranges::sort(state.navLayers, NavigationTarget::pathLengthOrder);

    // Print layer information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navLayers.size();
      os << " layer candidates found at path(s): ";
      for (auto& lc : state.navLayers) {
        os << lc.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navLayers.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No layer candidates found.");
    }
  }

  /// @brief Resolve compatible boundaries
  ///
  /// This function resolves the compatible boundaries for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveBoundaries(State& state, const Vector3& position,
                         const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible boundaries.");

    NavigationOptions<Surface> navOpts;
    navOpts.startObject = state.currentSurface;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    ACTS_VERBOSE(volInfo(state)
                 << "Try to find boundaries, we are at: " << toString(position)
                 << ", dir: " << toString(direction));

    if (m_geometryVersion == GeometryVersion::Gen1) {
      // Request the compatible boundaries
      state.navBoundaries = state.currentVolume->compatibleBoundaries(
          state.options.geoContext, position, direction, navOpts, logger());
      std::ranges::sort(state.navBoundaries, NavigationTarget::pathLengthOrder);
    } else {
      // Gen 3 !
      state.stream.reset();
      AppendOnlyNavigationStream appendOnly{state.stream};
      NavigationArguments args;
      args.position = position;
      args.direction = direction;
      args.wantsPortals = true;
      args.wantsSurfaces = false;
      state.currentVolume->initializeNavigationCandidates(args, appendOnly,
                                                          logger());

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << state.stream.candidates().size()
                   << " navigation candidates.");

      state.stream.initialize(state.options.geoContext, {position, direction},
                              BoundaryTolerance::None(),
                              state.options.surfaceTolerance);

      state.navBoundaries.clear();
      for (const NavigationTarget& candidate : state.stream.candidates()) {
        if (!detail::checkPathLength(candidate.pathLength(),
                                     state.options.nearLimit,
                                     state.options.farLimit, logger())) {
          continue;
        }

        state.navBoundaries.push_back(candidate);
      }
    }

    // Print boundary information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navBoundaries.size();
      os << " boundary candidates found at path(s): ";
      for (const auto& bc : state.navBoundaries) {
        os << bc.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navBoundaries.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No boundary candidates found.");
    }
  }

  /// @brief Check if the navigator is inactive
  ///
  /// This function checks if the navigator is inactive.
  ///
  /// @param state The navigation state
  ///
  /// @return True if the navigator is inactive
  bool inactive(const State& state) const {
    // Turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    if (state.navigationBreak) {
      return true;
    }

    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  // Cached so we don't have to query the TrackingGeometry constantly.
  TrackingGeometry::GeometryVersion m_geometryVersion{};

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
