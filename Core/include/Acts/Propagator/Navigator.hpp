// This file is part of the Acts project.
//
// Copyright (C) 2016-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <sstream>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief struct for the Navigation options that are forwarded to
///        the geometry
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

  /// Force intersection with boundaries
  bool forceIntersectBoundaries = false;
};

/// @brief Steers the propagation through the geometry by adjusting the step
///        size and providing the next surface to be targeted.
///
/// The Navigator is part of the propagation and responsible for steering
/// the step size in order to encounter all the relevant surfaces which are
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
  using Surfaces = std::vector<const Surface*>;

  using NavigationSurfaces =
      boost::container::small_vector<SurfaceIntersection, 10>;

  using NavigationLayers =
      boost::container::small_vector<LayerIntersection, 10>;

  using NavigationBoundaries =
      boost::container::small_vector<BoundaryIntersection, 4>;

  using ExternalSurfaces = std::multimap<std::uint64_t, GeometryIdentifier>;

  /// The navigation stage
  enum struct Stage : int {
    undefined = 0,
    surfaceTarget = 1,
    layerTarget = 2,
    boundaryTarget = 3
  };

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

  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    void insertExternalSurface(GeometryIdentifier geoid) {
      externalSurfaces.insert(
          std::pair<std::uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
    }

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface index of the navigation state
    std::size_t navSurfaceIndex = navSurfaces.size();

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer index of the navigation state
    std::size_t navLayerIndex = navLayers.size();

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary index of the navigation state
    std::size_t navBoundaryIndex = navBoundaries.size();

    SurfaceIntersection& navSurface() {
      return navSurfaces.at(navSurfaceIndex);
    }
    LayerIntersection& navLayer() { return navLayers.at(navLayerIndex); }
    BoundaryIntersection& navBoundary() {
      return navBoundaries.at(navBoundaryIndex);
    }

    /// Navigation state: the world volume
    const TrackingVolume* worldVolume = nullptr;

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the current layer
    const Layer* currentLayer = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;
    /// Navigation state: the target layer
    const Layer* targetLayer = nullptr;
    /// Navigation state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator if the target is reached
    bool targetReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;
    // The navigation stage (@todo: integrate break, target)
    Stage navigationStage = Stage::undefined;
    /// Force intersection with boundaries
    bool forceIntersectBoundaries = false;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {}

  State makeState(const Options& options) const {
    assert(options.startSurface != nullptr && "Start surface must be set");

    State state(options);
    state.startSurface = options.startSurface;
    state.targetSurface = options.targetSurface;
    return state;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool targetReached(const State& state) const { return state.targetReached; }

  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  void currentSurface(State& state, const Surface* surface) const {
    state.currentSurface = surface;
  }

  void targetReached(State& state, bool targetReached) const {
    state.targetReached = targetReached;
  }

  void navigationBreak(State& state, bool navigationBreak) const {
    state.navigationBreak = navigationBreak;
  }

  /// @brief Initialize call - start of navigation
  ///
  /// @param [in,out] state the navigation state
  void initialize(State& state, const Vector3& position,
                  const Vector3& direction) const {
    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    // Set the world volume if it is not set
    if (state.worldVolume == nullptr) {
      state.worldVolume = m_cfg.trackingGeometry->highestTrackingVolume();
    }

    // We set the current surface to the start surface
    // for eventual post-update action, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    state.currentSurface = state.startSurface;
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Current surface set to start surface "
                                  << state.currentSurface->geometryId());
    }

    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Surface.");
      // assign the current layer and volume by association
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
      // current volume and layer search through global search
      ACTS_VERBOSE(volInfo(state)
                   << "Starting from position " << toString(position)
                   << " and direction " << toString(direction));
      state.startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);
      state.startLayer = state.startVolume != nullptr
                             ? state.startVolume->associatedLayer(
                                   state.options.geoContext, position)
                             : nullptr;
      if (state.startVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      }
    }

    // Set the start volume as current volume
    state.currentVolume = state.startVolume;
    // Set the start layer as current layer
    state.currentLayer = state.startLayer;

    if (state.startLayer != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start layer to be resolved.");
      // We provide the layer to the resolve surface method in this case
      bool startResolved = resolveSurfaces(state, position, direction);
      if (!startResolved && state.startLayer == state.targetLayer) {
        ACTS_VERBOSE(volInfo(state)
                     << "Start is target layer and we have no surface "
                        "candidates. Nothing left to do.");
        // set the navigation break
        state.navigationBreak = true;
      }
    }
  }

  /// @brief Navigator estimateNextTarget call
  ///
  /// Call options
  /// (a) there are still surfaces to be resolved: handle those
  /// (b) there no surfaces but still layers to be resolved, handle those
  /// (c) there are no surfaces nor layers to be resolved, handle boundary
  ///
  /// @param [in,out] state the navigation state
  /// @param [in] position the current position
  /// @param [in] direction the current direction
  ///
  /// @return the next target surface intersection
  SurfaceIntersection estimateNextTarget(State& state, const Vector3& position,
                                         const Vector3& direction) const {
    // Check if the navigator is inactive
    if (inactive(state)) {
      return SurfaceIntersection::invalid();
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::estimateNextTarget.");

    // Reset the current surface
    state.currentSurface = nullptr;

    // Try targeting the surfaces - then layers - then boundaries
    if (state.navigationStage <= Stage::surfaceTarget &&
        targetSurfaces(state)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
      return state.navSurface();
    }

    if (state.navigationStage <= Stage::layerTarget &&
        targetLayers(state, position, direction)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
      return state.navLayer().first;
    }

    if (targetBoundaries(state, position, direction)) {
      ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
      return state.navBoundary().first;
    }

    ACTS_VERBOSE(volInfo(state)
                 << "No further navigation action, proceed to target.");
    // Set navigation break and release the navigation step size
    state.navigationBreak = true;
    return SurfaceIntersection::invalid();
  }

  void registerSurfaceStatus(State& state, const Vector3& position,
                             const Vector3& direction, const Surface& surface,
                             IntersectionStatus surfaceStatus) const {
    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    ACTS_VERBOSE(volInfo(state)
                 << "Entering Navigator::registerSurfaceStatus.");

    // Reset the current surface
    state.currentSurface = nullptr;

    if (surfaceStatus == IntersectionStatus::reachable) {
      ACTS_VERBOSE(volInfo(state) << "Stay focussed on surface.");
      return;
    }

    if (surfaceStatus == IntersectionStatus::onSurface) {
      state.currentSurface = &surface;
    }

    if (state.navigationStage <= Stage::surfaceTarget &&
        state.navSurfaceIndex < state.navSurfaces.size() &&
        state.navSurface().object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling surface status.");
      if (state.navSurfaceIndex < state.navSurfaces.size()) {
        // Set the navigation stage back to surface handling
        state.navigationStage = Stage::surfaceTarget;
        ++state.navSurfaceIndex;
      } else {
        // This was the last surface, switch to layers
        ACTS_VERBOSE(volInfo(state) << "Target layers.");
        state.navigationStage = Stage::layerTarget;
        ++state.navLayerIndex;
      }
      return;
    }

    if (state.navigationStage <= Stage::layerTarget &&
        state.navLayerIndex < state.navLayers.size() &&
        state.navLayer().first.object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling layer status.");

      if (surfaceStatus == IntersectionStatus::onSurface) {
        if (state.navLayerIndex < state.navLayers.size()) {
          state.currentLayer = state.navLayer().second;
          if (resolveSurfaces(state, position, direction)) {
            // Set the navigation stage back to surface handling
            state.navigationStage = Stage::surfaceTarget;
          } else {
            // Set the navigation stage to layer target
            state.navigationStage = Stage::layerTarget;
          }
          ++state.navLayerIndex;
        } else {
          // This was the last layer, switch to boundaries
          ACTS_VERBOSE(volInfo(state) << "Target boundaries.");
          state.navigationStage = Stage::boundaryTarget;
        }
      } else if (surfaceStatus == IntersectionStatus::missed) {
        // Set the navigation stage back to layer target
        ACTS_VERBOSE(volInfo(state) << "Staying focussed on layer.");
        state.navigationStage = Stage::layerTarget;
        ++state.navLayerIndex;
      }
      return;
    }

    if (state.navigationStage <= Stage::boundaryTarget &&
        state.navBoundaryIndex < state.navBoundaries.size() &&
        state.navBoundary().first.object() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling boundary status.");

      if (surfaceStatus == IntersectionStatus::onSurface) {
        // Update volume information: get the attached volume information
        const BoundarySurface* boundary = state.navBoundary().second;
        state.currentVolume = boundary->attachedVolume(state.options.geoContext,
                                                       position, direction);
        // reset
        state.navSurfaces.clear();
        state.navSurfaceIndex = state.navSurfaces.size();
        state.navLayers.clear();
        state.navLayerIndex = state.navLayers.size();
        state.navBoundaries.clear();
        state.navBoundaryIndex = state.navBoundaries.size();
        state.currentLayer = nullptr;
        if (state.currentVolume == nullptr) {
          // No volume anymore: end of known world
          ACTS_VERBOSE(
              volInfo(state)
              << "No more volume to progress to, stopping navigation.");
          // Navigation break & release navigation stepping
          state.navigationBreak = true;
        } else {
          ACTS_VERBOSE(volInfo(state) << "Volume updated.");
          state.navigationStage = Stage::undefined;
        }
      } else if (surfaceStatus == IntersectionStatus::missed) {
        // Set the navigation stage back to boundary target
        ACTS_VERBOSE(volInfo(state) << "Staying focussed on boundary.");
        state.navigationStage = Stage::boundaryTarget;
        ++state.navBoundaryIndex;
      }
      return;
    }
  }

 private:
  const SurfaceIntersection& candidateIntersection(
      const NavigationSurfaces& surfaces, std::size_t index) const {
    return surfaces.at(index);
  }
  const SurfaceIntersection& candidateIntersection(
      const NavigationLayers& surfaces, std::size_t index) const {
    return surfaces.at(index).first;
  }
  const SurfaceIntersection& candidateIntersection(
      const NavigationBoundaries& surfaces, std::size_t index) const {
    return surfaces.at(index).first;
  }

  /// Loop over surface candidates here:
  ///  - if an intersect is  valid but not yet reached
  ///    then return with updated step size
  ///  - if an intersect is not valid, switch to next
  ///
  /// @param [in,out] state the navigation state
  ///
  /// boolean return triggers exit to stepper
  bool targetSurfaces(State& state) const {
    if (state.navigationBreak) {
      return false;
    }

    // The call that we are on a layer and have not yet resolved the surfaces
    // No surfaces, do not return to stepper
    if (state.navSurfaces.empty() ||
        state.navSurfaceIndex == state.navSurfaces.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No surfaces present, target at layer first.");
      return false;
    }

    if (state.navSurfaceIndex < state.navSurfaces.size()) {
      return true;
    }

    // Reached the end of the surface iteration

    // first clear the surface cache
    state.navSurfaces.clear();
    state.navSurfaceIndex = state.navSurfaces.size();

    if (state.navLayerIndex != state.navLayers.size()) {
      ACTS_VERBOSE(volInfo(state)
                   << "Last surface on layer reached, switching layer.");
      // now switch to the next layer
      ++state.navLayerIndex;
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Last surface on layer reached, and no layer.");
      state.navigationBreak = (state.currentVolume == state.targetVolume);
    }

    // Do not return to the propagator
    return false;
  }

  /// @brief Target layer candidates.
  ///
  /// We are now trying to advance to the next layer (with surfaces)
  /// Check if we are on the representing surface of the layer pointed
  /// at by navLayerIndex. If so, we unpack the compatible surfaces
  /// (determined by straight line intersect), and set up the index
  /// so that the next postStep() call will enter the surface
  /// check mode above. If no surfaces are found, we skip the layer.
  /// If we unpack a surface, the step size is set to the path length
  /// to the first surface, as determined by straight line intersect.
  ///
  /// @param [in,out] state the navigation state
  ///
  /// @return boolean return triggers exit to stepper
  bool targetLayers(State& state, const Vector3& position,
                    const Vector3& direction) const {
    if (state.navigationBreak) {
      return false;
    }

    // if there are no layers, go back to the navigator (not stepper yet)
    if (state.navLayers.empty()) {
      ACTS_VERBOSE(volInfo(state)
                   << "No layers present, resolve volume first.");

      if (resolveLayers(state, position, direction)) {
        // The layer resolving worked
        return true;
      }
    }

    if (state.navLayerIndex < state.navLayers.size()) {
      return true;
    }

    // Screen output
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << "Last layer";
      if (state.currentVolume == state.targetVolume) {
        os << " (final volume) done, proceed to target.";
      } else {
        os << " done, target volume boundary.";
      }
      logger().log(Logging::VERBOSE, os.str());
    }
    // Set the navigation break if necessary
    state.navigationBreak = (state.currentVolume == state.targetVolume);
    return false;
  }

  /// @brief Navigation through volumes
  ///
  /// This is the boundary check routine. If the code above set up the
  /// boundary surface index, we advance through them here. If we are on
  /// the boundary surface, we set the current surface to the boundary
  /// surface, and get the volume pointed at by the boundary surface.  Next
  /// we unpack the layers from that volume. If the volume contains layers
  /// we set the step size to the straight line path length to the first
  /// layer.  If we don't find a next volume, the navigationBreak
  /// indicator is set.  This ends the navigation. Finally, the boundary
  /// index is cleared, so that the subsequent call goes back to
  /// the layer iteration logic.
  ///
  /// If we are not on the current boundary surface, we try the next one.
  /// The index is advanced and the step size is set. If no straight
  /// line intersect is found, the boundary surface is skipped.
  /// If we are out of boundary surfaces, the navigation is terminated.
  ///
  /// @param [in,out] state the navigation state
  ///
  /// boolean return triggers exit to stepper
  bool targetBoundaries(State& state, const Vector3& position,
                        const Vector3& direction) const {
    if (state.navigationBreak) {
      return false;
    }

    if (state.currentVolume == nullptr) {
      ACTS_VERBOSE(volInfo(state)
                   << "No sufficient information to resolve boundary, "
                      "stopping navigation.");
      return false;
    }

    if (state.currentVolume == state.targetVolume) {
      ACTS_VERBOSE(volInfo(state)
                   << "In target volume: no need to resolve boundary, "
                      "stopping navigation.");
      state.navigationBreak = true;
      return true;
    }

    // Helper function to find boundaries
    auto findBoundaries = [&]() -> void {
      // The navigation options
      NavigationOptions<Surface> navOpts;
      // Exclude the current surface in case it's a boundary
      navOpts.startObject = state.currentSurface;
      navOpts.nearLimit = state.options.nearLimit;
      navOpts.farLimit = state.options.farLimit;
      navOpts.forceIntersectBoundaries = state.forceIntersectBoundaries;

      ACTS_VERBOSE(volInfo(state)
                   << "Try to find boundaries, we are at: "
                   << toString(position) << ", dir: " << toString(direction));

      // Evaluate the boundary surfaces
      state.navBoundaries = state.currentVolume->compatibleBoundaries(
          state.options.geoContext, position, direction, navOpts, logger());
      std::sort(state.navBoundaries.begin(), state.navBoundaries.end(),
                [](const auto& a, const auto& b) {
                  return SurfaceIntersection::pathLengthOrder(a.first, b.first);
                });

      // Print boundary information
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << state.navBoundaries.size();
        os << " boundary candidates found at path(s): ";
        for (auto& bc : state.navBoundaries) {
          os << bc.first.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }

      // Set the begin index
      state.navBoundaryIndex = 0;
    };

    // No boundaries are assigned yet, find them
    if (state.navBoundaries.empty()) {
      findBoundaries();
    }

    if (state.navBoundaryIndex < state.navBoundaries.size()) {
      return true;
    }

    // We have to leave the volume somehow, so try again
    state.navBoundaries.clear();
    ACTS_VERBOSE(volInfo(state) << "Boundary navigation lost, re-targetting.");
    state.forceIntersectBoundaries = true;
    findBoundaries();
    if (!state.navBoundaries.empty()) {
      // Resetting intersection check for boundary surfaces
      state.forceIntersectBoundaries = false;
      return true;
    }

    // Tried our best, but couldn't do anything
    return false;
  }

  /// @brief Resolve the surfaces of this layer
  ///
  /// @param [in,out] state the navigation state
  ///
  /// boolean return triggers exit to stepper
  bool resolveSurfaces(State& state, const Vector3& position,
                       const Vector3& direction) const {
    // get the layer and layer surface
    const Layer* currentLayer = state.currentLayer;
    assert(currentLayer != nullptr && "Current layer is not set.");
    const Surface* layerSurface = &currentLayer->surfaceRepresentation();
    // Use navigation parameters and NavigationOptions
    NavigationOptions<Surface> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.currentSurface;
    navOpts.endObject = state.targetSurface;

    std::vector<GeometryIdentifier> externalSurfaces;
    if (!state.options.externalSurfaces.empty()) {
      auto layerID = layerSurface->geometryId().layer();
      auto externalSurfaceRange =
          state.options.externalSurfaces.equal_range(layerID);
      navOpts.externalSurfaces.reserve(
          state.options.externalSurfaces.count(layerID));
      for (auto itSurface = externalSurfaceRange.first;
           itSurface != externalSurfaceRange.second; itSurface++) {
        navOpts.externalSurfaces.push_back(itSurface->second);
      }
    }
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    // get the surfaces
    state.navSurfaces = currentLayer->compatibleSurfaces(
        state.options.geoContext, position, direction, navOpts);
    std::sort(state.navSurfaces.begin(), state.navSurfaces.end(),
              SurfaceIntersection::pathLengthOrder);

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

    // Set the index to the first
    state.navSurfaceIndex = 0;

    if (!state.navSurfaces.empty()) {
      // Surface candidates have been found
      return true;
    }

    ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
    return false;
  }

  /// @brief Navigation through layers
  ///
  /// Resolve layers.
  ///
  /// This initializes the layer candidates when starting
  /// or when entering a new volume
  ///
  /// @param [in,out] state the navigation state
  ///
  /// @return boolean return triggers exit to stepper
  bool resolveLayers(State& state, const Vector3& position,
                     const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

    // Create the navigation options
    // - and get the compatible layers, start layer will be excluded
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
    std::sort(state.navLayers.begin(), state.navLayers.end(),
              [](const auto& a, const auto& b) {
                return SurfaceIntersection::pathLengthOrder(a.first, b.first);
              });

    // Print layer information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navLayers.size();
      os << " layer candidates found at path(s): ";
      for (auto& lc : state.navLayers) {
        os << lc.first.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    // Set the index to the first
    state.navLayerIndex = 0;

    if (!state.navLayers.empty()) {
      // Layer candidates have been found
      return true;
    }

    // Screen output - no layer candidates found
    ACTS_VERBOSE(volInfo(state) << "No compatible layer candidates found.");
    return false;
  }

  /// Inactive
  ///
  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @param [in,out] state the navigation state
  ///
  /// boolean return triggers exit to stepper
  bool inactive(const State& state) const {
    // Void behavior in case no tracking geometry is present
    if (m_cfg.trackingGeometry == nullptr) {
      return true;
    }

    // turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    // Navigation break handling
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

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
