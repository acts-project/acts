// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/AnyIntersection.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <cstddef>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief struct for the Navigation options that are forwarded to
///        the geometry
///
/// @tparam propagator_state_t Type of the object for navigation state
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryCheck boundaryCheck = BoundaryCheck(true);

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// object to check against: at start
  const object_t* startObject = nullptr;
  /// object to check against: at end
  const object_t* endObject = nullptr;

  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The minimum distance for a surface to be considered
  double nearLimit = 0;
  /// The maximum distance for a surface to be considered
  double farLimit = std::numeric_limits<double>::max();
};

/// Navigator class
///
/// This is an Actor to be added to the ActorList in order to navigate
/// through the static tracking geometry setup.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation  flow is hit, it is
/// set to the propagation satate, such that other actors can deal with it.
/// This navigation actor thus always needs to run first!
/// It does two things: it figures out the order of volumes, layers and
/// surfaces.
///
/// The current target surface is the surface pointed to by of the index
/// for the surfaces, layers or volume boundaries.
/// If a surface is found, the state.navigation.currentSurface
/// pointer is set. This  enables subsequent actors to react. Secondly, this
/// actor uses the ordered indices to figure out which surface, layer or
/// volume boundary is _supposed_ to be hit next. It then sets the maximum
/// step size to the path length found out by straight line intersection. If
/// the state is not on surface, it also  re-computes the step size, to make
/// sure we end up at the desired surface.
///
class Navigator {
 public:
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// Configuration for this Navigator
    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;

    /// Whether to perform boundary checks for layer resolving (improves
    /// navigation for bended tracks)
    BoundaryCheck boundaryCheckLayerResolving = BoundaryCheck(true);
  };

  struct NavigationCandidate {
    AnyIntersection intersection;
    BoundaryCheck boundaryCheck = BoundaryCheck(true);

    NavigationCandidate(AnyIntersection _intersection,
                        BoundaryCheck _boundaryCheck)
        : intersection(std::move(_intersection)),
          boundaryCheck(std::move(_boundaryCheck)) {}

    static bool forwardOrder(const NavigationCandidate& aCandidate,
                             const NavigationCandidate& bCandidate) {
      return Intersection3D::forwardOrder(
          aCandidate.intersection.intersection(),
          bCandidate.intersection.intersection());
    }
  };

  using Surfaces = std::vector<const Surface*>;

  using NavigationCandidates =
      boost::container::small_vector<NavigationCandidate, 24>;

  using ExternalSurfaces = std::multimap<uint64_t, GeometryIdentifier>;

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    /// the vector of navigation surfaces to work through
    NavigationCandidates candidates = {};
    /// the current surface index of the navigation state
    int candidateIndex = 0;

    auto candidate() const { return candidates.at(candidateIndex); }

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    /// Navigation state: the world volume
    const TrackingVolume* worldVolume = nullptr;

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state - the current layer
    const Layer* currentLayer = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target volume
    const TrackingVolume* targetVolume = nullptr;
    /// Navigation state: the target layer
    const Layer* targetLayer = nullptr;
    /// Navigation state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator for start layer treatment
    bool startLayerResolved = false;
    /// Indicator if the target is reached
    bool targetReached = false;
    /// Indicator that the last VolumeHierarchy surface was reached
    /// skip the next layer targeting to the next boundary/volume
    bool lastHierarchySurfaceReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;

    /// Navigation reset requested
    bool resetRequested = false;
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {}

  State makeState(const Surface* startSurface,
                  const Surface* targetSurface) const {
    State result;
    result.startSurface = startSurface;
    result.targetSurface = targetSurface;
    return result;
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

  void insertExternalSurface(State& state, GeometryIdentifier geoid) const {
    state.externalSurfaces.insert(
        std::pair<uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
  }

  /// @brief Initialize call - start of propagation
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, const stepper_t& stepper) const {
    ACTS_VERBOSE("Initialization.");

    if (state.navigation.startSurface &&
        state.navigation.startSurface->associatedLayer()) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Surface.");
      state.navigation.startLayer =
          state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume =
          state.navigation.startLayer->trackingVolume();
    } else if (state.navigation.startVolume) {
      ACTS_VERBOSE(
          "Fast start initialization through association from Volume.");
      state.navigation.startLayer =
          state.navigation.startVolume->associatedLayer(
              state.geoContext, stepper.position(state.stepping));
    } else {
      ACTS_VERBOSE("Slow start initialization through search.");
      ACTS_VERBOSE("Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
      state.navigation.startVolume =
          m_cfg.trackingGeometry->lowestTrackingVolume(
              state.geoContext, stepper.position(state.stepping));
      state.navigation.startLayer =
          state.navigation.startVolume
              ? state.navigation.startVolume->associatedLayer(
                    state.geoContext, stepper.position(state.stepping))
              : nullptr;
    }

    // Special handling if the start surface is a boundary surface since the
    // start volume might not be correct.
    if (state.navigation.startSurface != nullptr) {
      for (const auto& boundary :
           state.navigation.startVolume->boundarySurfaces()) {
        if (state.navigation.startSurface ==
            &boundary->surfaceRepresentation()) {
          state.navigation.startVolume = boundary->attachedVolume(
              state.geoContext, stepper.position(state.stepping),
              state.options.direction * stepper.direction(state.stepping));
          assert(state.navigation.startVolume != nullptr &&
                 "Start volume not resolved.");
          ACTS_VERBOSE("We are starting from a boundary surface "
                       << state.navigation.startSurface->geometryId()
                       << ". Attached volume is "
                       << state.navigation.startVolume->geometryId());
          break;
        }
      }
    }

    // Initialize current volume, layer and surface
    {
      state.navigation.worldVolume =
          m_cfg.trackingGeometry->highestTrackingVolume();

      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.currentVolume) {
        ACTS_VERBOSE(volInfo(state) << "Start volume resolved.");
      } else {
        ACTS_ERROR("Start volume not resolved.");
      }

      state.navigation.currentLayer = state.navigation.startLayer;
      if (state.navigation.currentLayer) {
        ACTS_VERBOSE(volInfo(state)
                     << "Start layer resolved to "
                     << state.navigation.currentLayer->geometryId() << " .");
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start layer set.");
      }

      state.navigation.currentSurface = state.navigation.startSurface;
      if (state.navigation.currentSurface) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current surface set to start surface "
                     << state.navigation.currentSurface->geometryId());
      } else {
        ACTS_VERBOSE(volInfo(state) << "No start surface set.");
      }
    }

    // Initialize navigation candidates for the start volume
    {
      state.navigation.candidates.clear();
      state.navigation.candidateIndex = 0;

      initializeVolume(state, stepper);

      state.navigation.currentLayer = state.navigation.startLayer;
      if (state.navigation.currentLayer) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current layer set to start layer "
                     << state.navigation.currentLayer->geometryId());
      }

      initializeLayer(state, stepper);

      std::sort(state.navigation.candidates.begin(),
                state.navigation.candidates.end(),
                NavigationCandidate::forwardOrder);
    }
  }

  /// @brief Navigator pre step call
  ///
  /// Call options
  /// (a) there are still surfaces to be resolved: handle those
  /// (b) there no surfaces but still layers to be resolved, handle those
  /// (c) there are no surfaces nor layers to be resolved, handle boundary
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void preStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Check for reset
    if (state.navigation.resetRequested) {
      initialize(state, stepper);
      state.navigation.resetRequested = false;
    }

    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::preStep.");

    // Loop over the remaining navigation candidates
    for (; state.navigation.candidateIndex !=
           (int)state.navigation.candidates.size();
         ++state.navigation.candidateIndex) {
      // Screen output how much is left to try
      ACTS_VERBOSE(volInfo(state)
                   << (state.navigation.candidates.size() -
                       state.navigation.candidateIndex)
                   << " out of " << state.navigation.candidates.size()
                   << " surfaces remain to try.");

      const auto& candidate = state.navigation.candidate();
      const auto& intersection = candidate.intersection;
      const Surface* surface = intersection.representation();

      if (intersection.template checkType<SurfaceIntersection>()) {
        if (surface->associatedLayer() != state.navigation.currentLayer) {
          ACTS_VERBOSE(volInfo(state) << "Skip surface from a different layer");
          continue;
        }
      }

      // Screen output which surface you are on
      ACTS_VERBOSE(volInfo(state) << "Next surface candidate will be "
                                  << surface->geometryId());
      // Get the boundary check
      BoundaryCheck boundaryCheck = candidate.boundaryCheck;
      // Estimate the surface status
      auto surfaceStatus = stepper.updateSurfaceStatus(
          state.stepping, *surface, intersection.index(),
          state.options.direction, boundaryCheck,
          state.options.surfaceTolerance, logger());
      if (surfaceStatus == Intersection3D::Status::reachable) {
        ACTS_VERBOSE(volInfo(state)
                     << "Surface reachable, step size updated to "
                     << stepper.outputStepSize(state.stepping));
        break;
      }
    }

    if (state.navigation.candidateIndex ==
        (int)state.navigation.candidates.size()) {
      ACTS_ERROR(volInfo(state) << "Exhausted navigation candidates.");
      // Set navigation break and release the navigation step size
      state.navigation.navigationBreak = true;
      stepper.releaseStepSize(state.stepping);
    }

    // Navigator target always resets the current surface
    state.navigation.currentSurface = nullptr;
  }

  /// @brief Navigator post step call, will be called in two modes
  ///
  /// (a) It initializes the Navigation stream if start volume is
  ///     not yet defined:
  ///  - initialize the volume
  ///  - establish the start layer and start volume
  ///  - set the current surface to the start surface
  ///
  /// (b) It establishes the currentSurface status during
  ///     the propagation flow, currentSurface can be
  ///  - surfaces still to be handled within a layer
  ///  - layers still to be handled within a volume
  ///  - boundaries still to be handled to exit a volume
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void postStep(propagator_state_t& state, const stepper_t& stepper) const {
    // Check if the navigator is inactive
    if (inactive(state)) {
      return;
    }

    // Call the navigation helper prior to actual navigation
    ACTS_VERBOSE(volInfo(state) << "Entering navigator::postStep.");

    assert(state.navigation.currentSurface == nullptr &&
           "Current surface must be reset.");
    surfaceStatus(state, stepper);

    if (state.navigation.currentSurface == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Staying focussed on surface.");
      return;
    }

    ACTS_VERBOSE(volInfo(state)
                 << "Handle surface "
                 << state.navigation.currentSurface->geometryId());

    const auto& intersection = state.navigation.candidate().intersection;

    if (intersection.template checkType<SurfaceIntersection>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a surface");
    } else if (intersection.template checkType<LayerIntersection>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a layer");

      state.navigation.currentLayer =
          intersection.template object<LayerIntersection>();

      initializeLayer(state, stepper);

      std::sort(state.navigation.candidates.begin() +
                    state.navigation.candidateIndex + 1,
                state.navigation.candidates.end(),
                NavigationCandidate::forwardOrder);
    } else if (intersection.template checkType<BoundaryIntersection>()) {
      ACTS_VERBOSE(volInfo(state) << "This is a boundary");

      const auto* boundary =
          intersection.template object<BoundaryIntersection>();

      state.navigation.currentVolume = boundary->attachedVolume(
          state.geoContext, stepper.position(state.stepping),
          state.options.direction * stepper.direction(state.stepping));

      ACTS_VERBOSE(volInfo(state) << "Switched volume");

      state.navigation.candidates.clear();
      state.navigation.candidateIndex = 0;

      initializeVolume(state, stepper);

      std::sort(state.navigation.candidates.begin(),
                state.navigation.candidates.end(),
                NavigationCandidate::forwardOrder);
    } else {
      ACTS_ERROR(volInfo(state) << "Unknown intersection type");
    }
  }

 private:
  /// @brief Status call for test surfaces (surfaces, layers, boundaries)
  ///
  /// If there are surfaces to be handled, check if the current
  /// state is on the surface
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  /// @tparam navigation_surfaces_t Type of the propagator
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void surfaceStatus(propagator_state_t& state,
                     const stepper_t& stepper) const {
    const auto& candidate = state.navigation.candidate();
    const auto& intersection = candidate.intersection;
    // Take the current surface
    const Surface* surface = intersection.representation();
    // Check if we are at a surface
    // If we are on the surface pointed at by the index, we can make
    // it the current one to pass it to the other actors
    auto surfaceStatus = stepper.updateSurfaceStatus(
        state.stepping, *surface, intersection.index(), state.options.direction,
        candidate.boundaryCheck, state.options.surfaceTolerance, logger());
    if (surfaceStatus == Intersection3D::Status::onSurface) {
      ACTS_VERBOSE(volInfo(state)
                   << "Status Surface successfully hit, storing it.");
      // Set in navigation state, so actors and aborters can access it
      state.navigation.currentSurface = surface;
      if (state.navigation.currentSurface) {
        ACTS_VERBOSE(volInfo(state)
                     << "Current surface set to surface "
                     << state.navigation.currentSurface->geometryId());
      }
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  void initializeVolume(propagator_state_t& state,
                        const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Initialize volume");

    if (state.navigation.currentVolume == nullptr) {
      state.navigation.navigationBreak = true;
      ACTS_VERBOSE(volInfo(state) << "No volume set. Good luck.");
      return;
    }

    // Find boundary candidates
    {
      ACTS_VERBOSE(volInfo(state) << "Searching for compatible boundaries.");

      // The navigation options
      NavigationOptions<Surface> navOpts;
      // Exclude the current surface in case it's a boundary
      navOpts.startObject = state.navigation.currentSurface;
      navOpts.nearLimit = state.options.surfaceTolerance;
      navOpts.farLimit = std::numeric_limits<double>::max();

      ACTS_VERBOSE(volInfo(state)
                   << "Try to find boundaries, we are at: "
                   << stepper.position(state.stepping).transpose() << ", dir: "
                   << stepper.direction(state.stepping).transpose());

      auto boundaries = state.navigation.currentVolume->compatibleBoundaries(
          state.geoContext, stepper.position(state.stepping),
          state.options.direction * stepper.direction(state.stepping), navOpts,
          logger());

      // Screen output where they are
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << boundaries.size();
        os << " boundary candidates found at path(s): ";
        for (const auto& boundary : boundaries) {
          os << boundary.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }

      for (const auto& boundary : boundaries) {
        state.navigation.candidates.emplace_back(AnyIntersection(boundary),
                                                 BoundaryCheck(true));
      }
    }

    // Find layer candidates
    {
      ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

      // Create the navigation options
      // - and get the compatible layers, start layer will be excluded
      NavigationOptions<Layer> navOpts;
      navOpts.boundaryCheck = m_cfg.boundaryCheckLayerResolving;
      navOpts.resolveSensitive = m_cfg.resolveSensitive;
      navOpts.resolveMaterial = m_cfg.resolveMaterial;
      navOpts.resolvePassive = m_cfg.resolvePassive;
      navOpts.startObject =
          (state.navigation.currentVolume == state.navigation.startVolume)
              ? state.navigation.startLayer
              : nullptr;
      navOpts.nearLimit = state.options.surfaceTolerance;
      navOpts.farLimit = std::numeric_limits<double>::max();

      // Request the compatible layers
      auto layers = state.navigation.currentVolume->compatibleLayers(
          state.geoContext, stepper.position(state.stepping),
          state.options.direction * stepper.direction(state.stepping), navOpts);

      // Screen output where they are
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << layers.size();
        os << " layer candidates found at path(s): ";
        for (const auto& layer : layers) {
          os << layer.pathLength() << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }

      for (const auto& layer : layers) {
        state.navigation.candidates.emplace_back(
            AnyIntersection(layer), m_cfg.boundaryCheckLayerResolving);
      }
    }
  }

  template <typename propagator_state_t, typename stepper_t>
  auto initializeLayer(propagator_state_t& state,
                       const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Initialize layer");

    if (state.navigation.currentLayer == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "No layer set.");
      return;
    }

    auto surfaces = resolveSurfaces(state, stepper);

    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << surfaces.size();
      os << " surface candidates found at path(s): ";
      for (const auto& surface : surfaces) {
        os << surface.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    for (const auto& surface : surfaces) {
      state.navigation.candidates.emplace_back(AnyIntersection(surface),
                                               BoundaryCheck(true));
    }
  }

  /// @brief Resolve the surfaces of this layer
  ///
  /// @tparam propagator_state_t The state type of the propagator
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  auto resolveSurfaces(propagator_state_t& state,
                       const stepper_t& stepper) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible surfaces.");

    const Layer* layer = state.navigation.currentLayer;
    assert(layer != nullptr && "Current layer must be set.");

    // Use navigation parameters and NavigationOptions
    NavigationOptions<Surface> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = (layer == state.navigation.startLayer)
                              ? state.navigation.startSurface
                              : nullptr;
    navOpts.endObject = state.navigation.targetSurface;
    navOpts.nearLimit = state.options.surfaceTolerance;
    navOpts.farLimit = std::numeric_limits<double>::max();

    std::vector<GeometryIdentifier> externalSurfaces;
    if (!state.navigation.externalSurfaces.empty()) {
      auto layerID = layer->geometryId().layer();
      auto externalSurfaceRange =
          state.navigation.externalSurfaces.equal_range(layerID);
      navOpts.externalSurfaces.reserve(
          state.navigation.externalSurfaces.count(layerID));
      for (auto itSurface = externalSurfaceRange.first;
           itSurface != externalSurfaceRange.second; itSurface++) {
        navOpts.externalSurfaces.push_back(itSurface->second);
      }
    }

    // get the surfaces
    return layer->compatibleSurfaces(
        state.geoContext, stepper.position(state.stepping),
        state.options.direction * stepper.direction(state.stepping), navOpts);
  }

  /// Inactive
  ///
  /// This checks if a navigation break had been triggered or navigator
  /// is misconfigured
  ///
  /// @tparam propagator_state_t The state type of the propagator
  ///
  /// @param [in,out] state is the propagation state object
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t>
  bool inactive(propagator_state_t& state) const {
    // Void behavior in case no tracking geometry is present
    if (!m_cfg.trackingGeometry) {
      return true;
    }
    // turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }
    // handle navigation break
    if (state.navigation.navigationBreak) {
      return true;
    }
    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.navigation.currentVolume
                ? state.navigation.currentVolume->volumeName()
                : "No Volume") +
           " | ";
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
