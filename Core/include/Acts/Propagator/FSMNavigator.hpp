// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationOptions.hpp"
#include "Acts/Utilities/FiniteStateMachine.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string_view>
#include <tuple>

namespace Acts {

class FSMNavigator {
 public:
  struct states {
    struct Init {
      constexpr static std::string_view name = "Init";
    };

    struct VolumeToVolume {
      constexpr static std::string_view name = "VolumeToVolume";
    };

    struct LayerToLayer {
      constexpr static std::string_view name = "LayerToLayer";
    };

    struct SurfaceToSurface {
      constexpr static std::string_view name = "SurfaceToSurface";
    };

    struct ToBoundarySurface {
      constexpr static std::string_view name = "ToBoundarySurface";
    };

    struct Finished {
      constexpr static std::string_view name = "Finished";
    };
  };

  struct events {
    struct Status {
      constexpr static std::string_view name = "Status";
    };
    struct Target {
      constexpr static std::string_view name = "Target";
    };
    struct ResolveLayers {
      constexpr static std::string_view name = "ResolveLayers";
    };
    struct ResolveSurfaces {
      constexpr static std::string_view name = "ResolveSurfaces";
    };
    struct ResolveBoundarySurfaces {
      constexpr static std::string_view name = "ResolveBoundarySurfaces";
    };
  };

  using NavigationSurfaces = std::vector<SurfaceIntersection>;
  using NavigationSurfaceIter = NavigationSurfaces::iterator;

  using NavigationLayers = std::vector<LayerIntersection>;
  using NavigationLayerIter = NavigationLayers::iterator;

  using NavigationBoundaries = std::vector<BoundaryIntersection>;
  using NavigationBoundaryIter = NavigationBoundaries::iterator;

  using ExternalSurfaces = std::multimap<uint64_t, GeometryIdentifier>;

  struct State
      : FiniteStateMachine<State, states::Init, states::VolumeToVolume,
                           states::LayerToLayer, states::SurfaceToSurface,
                           states::ToBoundarySurface, states::Finished> {
    /// Externally provided surfaces - these are tried to be hit
    std::multimap<const Layer*, const Surface*> externalSurfaces = {};

    /// Navigation sate: the world volume
    const TrackingVolume* worldVolume = nullptr;

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
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
    // The navigation stage (@todo: integrate break, target)
    // Stage navigationStage = Stage::undefined;

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface iterator of the navigation state
    NavigationSurfaceIter navSurfaceIter = navSurfaces.end();

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer iterator of the navigation state
    NavigationLayerIter navLayerIter = navLayers.end();

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary iterator of the navigation state
    NavigationBoundaryIter navBoundaryIter = navBoundaries.end();

    friend FSMNavigator;
    friend fsm_base;

    // everything below is private to the FSMNavigator, and shouldn't be touched
    // from the outside
   private:
    std::string volstr(const TrackingVolume& tv) const {
      std::stringstream ss;
      ss << tv.volumeName() << " " << tv.geometryId();
      return ss.str();
    }

    // template <typename propagator_state_t, typename corrector_t>
    // void updateStep(propagator_state_t& state, const corrector_t&
    // [>navCorr<], double navigationStep, bool release = false) const {
    // //  update the step
    // state.stepping.stepSize.update(navigationStep,
    // detail::ConstrainedStep::actor, release);
    // ACTS_DEBUG("Navigation stepSize " << (release ? "released and " : "")
    // << "updated to "
    // << state.stepping.stepSize.toString());
    ///// If we have an initial step and are configured to modify it
    // if (state.stepping.pathAccumulated == 0.
    // and navCorr(state.stepping.stepSize)) {
    // debugLog(state, [&] {
    // std::stringstream dstream;
    // dstream << "Initial navigation step modified to ";
    // dstream << state.stepping.stepSize.toString();
    // return dstream.str();
    //});
    //}
    // }

    // Finite State event handlers and state enter/exit callbacks

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::Init&, const events::Status&,
                          const FSMNavigator& navigator,
                          propagator_state_t& state, const stepper_t& stepper) {
      if (!navigator.m_cfg.trackingGeometry) {
        return std::nullopt;
      }

      // get the world volume
      state.navigation.worldVolume =
          navigator.m_cfg.trackingGeometry->highestTrackingVolume();

      // in any case, we start from the start surface
      state.navigation.currentSurface = state.navigation.startSurface;
      if (state.navigation.currentSurface != nullptr) {
        ACTS_DEBUG("Current surface set to start surface: "
                   << state.navigation.currentSurface->geometryId());
      }

      // we need to figure out the current volume.
      // If we're on a layer, then we can ask the layer what it's volume is
      if (state.navigation.startSurface &&
          state.navigation.startSurface->associatedLayer()) {
        ACTS_DEBUG("Fast start initialization through layer association");

        state.navigation.startLayer =
            state.navigation.startSurface->associatedLayer();
        state.navigation.startVolume =
            state.navigation.startLayer->trackingVolume();
        // is also the current volume
        state.navigation.currentVolume = state.navigation.startVolume;
      } else {
        ACTS_DEBUG(
            "Slow start through volume hierarchy search for lowest volume");
        ACTS_DEBUG("Starting from position "
                   << toString(stepper.position(state.stepping))
                   << " and direction "
                   << toString(stepper.direction(state.stepping)));
        state.navigation.startVolume =
            navigator.m_cfg.trackingGeometry->lowestTrackingVolume(
                state.geoContext, stepper.position(state.stepping));
        state.navigation.startLayer =
            state.navigation.startVolume
                ? state.navigation.startVolume->associatedLayer(
                      state.geoContext, stepper.position(state.stepping))
                : nullptr;
        // Set the start volume as current volume
        state.navigation.currentVolume = state.navigation.startVolume;
        if (state.navigation.startVolume) {
          ACTS_DEBUG("Start volume resolved: "
                     << state.navigation.startVolume->volumeName() << " "
                     << state.navigation.startVolume->geometryId());
        } else {
          ACTS_ERROR("Unable to resolve start volume");
          return Terminated{};
        }
      }

      return std::nullopt;
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::Init&, const events::Target&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state,
                          const stepper_t& /*stepper*/) {
      ACTS_DEBUG("Initial target call.");

      // do we have a target surface?
      if (state.navigation.targetSurface != nullptr) {
        ACTS_DEBUG("Target surface is given: "
                   << state.navigation.targetSurface->geometryId());
      } else {
        ACTS_DEBUG("No target surface given");
      }

      return states::VolumeToVolume{};
    }

    template <typename propagator_state_t, typename stepper_t>
    void on_enter(const states::VolumeToVolume&, const FSMNavigator& navigator,
                  propagator_state_t& state, const stepper_t& stepper) {
      ACTS_DEBUG("Entering volume: "
                 << state.navigation.currentVolume->volumeName() << " "
                 << state.navigation.currentVolume->geometryId());

      stepper.releaseStepSize(state.stepping);

      // figure out the volume type, currently, we only support layer based
      // volumes
      dispatch(events::ResolveLayers{}, navigator, state, stepper);
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::VolumeToVolume&,
                          const events::ResolveLayers&,
                          const FSMNavigator& navigator,
                          propagator_state_t& state, const stepper_t& stepper) {
      ACTS_DEBUG("Resolving layers for volume "
                 << volstr(*state.navigation.currentVolume));

      // Check if we are in the start volume
      auto _startLayer =
          (state.navigation.currentVolume == state.navigation.startVolume)
              ? state.navigation.startLayer
              : nullptr;

      // Create the navigation options
      // - and get the compatible layers, start layer will be excluded
      NavigationOptions<Layer> navOpts(
          state.stepping.navDir, true, navigator.m_cfg.resolveSensitive,
          navigator.m_cfg.resolveMaterial, navigator.m_cfg.resolvePassive,
          _startLayer, nullptr);
      navOpts.pathLimit =
          state.stepping.stepSize.value(ConstrainedStep::aborter);

      state.navigation.navLayers =
          state.navigation.currentVolume->compatibleLayers(
              state.geoContext, stepper.position(state.stepping),
              stepper.direction(state.stepping), navOpts,
              LoggerWrapper{logger()});

      ACTS_DEBUG("Found " << state.navigation.navLayers.size()
                          << " layer candidates");

      if (state.navigation.navLayers.empty()) {
        return states::ToBoundarySurface{};
      }

      state.navigation.navLayerIter = state.navigation.navLayers.begin();

      if (logger().doPrint(Logging::DEBUG)) {
        std::stringstream ss;
        ss << "Layer candidates found at path(s): ";
        for (const auto& [ix, layer, surface, dir] :
             state.navigation.navLayers) {
          ss << ix.pathLength << " ";
        }
        ACTS_DEBUG(ss.str());
      }

      stepper.updateStepSize(state.stepping, *state.navigation.navLayerIter,
                             true);

      return states::LayerToLayer{};
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::LayerToLayer&, const events::Status&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      // did we hit the layer?
      const Surface* layerSurface =
          state.navigation.navLayerIter->representation;

      auto targetStatus =
          stepper.updateSurfaceStatus(state.stepping, *layerSurface, true);

      if (targetStatus == Intersection3D::Status::onSurface) {
        ACTS_DEBUG("Layer hit, current surface is now: "
                   << layerSurface->geometryId());
        state.navigation.currentSurface = layerSurface;
        return std::nullopt;  // stay in state
      } else {
        ACTS_DEBUG("Layer not hit. Continue");
        state.navigation.currentSurface = nullptr;
        return std::nullopt;  // stay in state
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::LayerToLayer&, const events::Target&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      // did we hit the layer?
      auto& [init_ix, layer, surface, dir] = *state.navigation.navLayerIter;
      if (surface == state.navigation.currentSurface) {
        return states::SurfaceToSurface{};
      } else {
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);

        while (state.navigation.navLayerIter !=
               state.navigation.navLayers.end()) {
          const Surface* layerSurface =
              state.navigation.navLayerIter->representation;

          auto layerStatus =
              stepper.updateSurfaceStatus(state.stepping, *layerSurface, true);
          if (layerStatus != Intersection3D::Status::reachable) {
            ACTS_DEBUG("During approach of layer "
                       << layer->geometryId()
                       << " intersection became unreachable => skipping layer "
                          "candidate");
            ++state.navigation.navLayerIter;
            continue;
          } else {
            // update straight line estimation
            ACTS_DEBUG("Proceeding towards layer: "
                       << layerSurface->geometryId() << ", updated step size: "
                       << stepper.outputStepSize(state.stepping));
            return std::nullopt;  // stay in state
          }
        }

        if (state.navigation.navLayerIter == state.navigation.navLayers.end()) {
          ACTS_DEBUG("No further layers in volume, target boundary surfaces");
          // state.navigation.navLayers.clear();
          // state.navigation.navLayerIter = stat.navigation.navLayers.end();
          return states::ToBoundarySurface{};
        }

        // we shouldn't get here
        return Terminated{};
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    void on_enter(const states::LayerToLayer&, const FSMNavigator& navigator,
                  propagator_state_t& state, const stepper_t& stepper) {
      ACTS_DEBUG("Entering layer to layer");

      if (state.navigation.navLayerIter == state.navigation.navLayers.end()) {
        ACTS_DEBUG("No further layers in volume, target boundary surfaces");
        // state.navigation.navLayers.clear();
        // state.navigation.navLayerIter = stat.navigation.navLayers.end();
        // @TODO: This needs to transition to ToBoundarySurface right away
        // dispatch(states::ToBoundarySurface{}, navigator, state, stepper);
        setState(states::ToBoundarySurface{}, navigator, state, stepper);
      } else {
        // make sure we're targeting that one
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);

        while (state.navigation.navLayerIter !=
               state.navigation.navLayers.end()) {
          auto& [init_ix, layer, layerSurface, dir] =
              *state.navigation.navLayerIter;

          auto layerStatus =
              stepper.updateSurfaceStatus(state.stepping, *layerSurface, true);
          if (layerStatus != Intersection3D::Status::reachable) {
            ACTS_DEBUG("During approach of layer "
                       << layer->geometryId()
                       << " intersection became unreachable => skipping layer "
                          "candidate");
            ++state.navigation.navLayerIter;
            continue;
          } else {
            // update straight line estimation
            ACTS_DEBUG("Proceeding towards layer: "
                       << layerSurface->geometryId() << ", updated step size: "
                       << stepper.outputStepSize(state.stepping));
            return;  // stay in state
          }
        }

        if (state.navigation.navLayerIter == state.navigation.navLayers.end()) {
          ACTS_DEBUG("No further layers in volume, target boundary surfaces");
          setState(states::ToBoundarySurface{}, navigator, state, stepper);
        }
      }

      // dispatch(events::ResolveSurfaces{}, navigator, state, stepper);
    }

    template <typename propagator_state_t, typename stepper_t>
    void on_enter(const states::SurfaceToSurface&,
                  const FSMNavigator& navigator, propagator_state_t& state,
                  const stepper_t& stepper) {
      ACTS_DEBUG("Entering surface to surface on layer "
                 << state.navigation.navLayerIter->object->geometryId());

      dispatch(events::ResolveSurfaces{}, navigator, state, stepper);
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::SurfaceToSurface&,
                          const events::ResolveSurfaces&,
                          const FSMNavigator& navigator,
                          propagator_state_t& state, const stepper_t& stepper) {
      const Layer* navLayer = state.navigation.navLayerIter->object;

      ACTS_DEBUG("Resolving surfaces for layer: " << navLayer->geometryId());

      const Surface* layerSurface = state.navigation.currentSurface;

      // are we on the start layer
      bool onStart = (navLayer == state.navigation.startLayer);

      NavigationOptions<Surface> navOpts(
          state.stepping.navDir, true, navigator.config().resolveSensitive,
          navigator.config().resolveMaterial, navigator.config().resolvePassive,
          layerSurface, state.navigation.targetSurface);

      // Check the limit
      navOpts.pathLimit =
          state.stepping.stepSize.value(ConstrainedStep::aborter);

      // No overstepping on start layer, otherwise ask the stepper
      navOpts.overstepLimit = onStart ? s_onSurfaceTolerance
                                      : stepper.overstepLimit(state.stepping);

      state.navigation.navSurfaces = navLayer->compatibleSurfaces(
          state.geoContext, stepper.position(state.stepping),
          stepper.direction(state.stepping), navOpts);

      stepper.releaseStepSize(state.stepping);

      // the number of layer candidates
      if (!state.navigation.navSurfaces.empty()) {
        if (logger().doPrint(Logging::VERBOSE)) {
          std::ostringstream os;
          os << state.navigation.navSurfaces.size();
          os << " surface candidates found at path(s): \n";
          for (auto& sfc : state.navigation.navSurfaces) {
            os << " - " << sfc.object->geometryId() << " at "
               << sfc.intersection.pathLength << "\n ";
          }
          logger().log(Logging::VERBOSE, os.str());
        }

        // set the iterator
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": targeting next surface: "
                     << state.navigation.navSurfaceIter->object->geometryId())
        // The stepper updates the step size ( single / multi component)
        stepper.updateStepSize(state.stepping, *state.navigation.navSurfaceIter,
                               true);
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": Navigation stepSize updated to "
                     << stepper.outputStepSize(state.stepping));
        return std::nullopt;  // stay in surface to surface

      } else {
        state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": No surface candidates found.");
        // state.navigation.navSurfaces.clear();
        // state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
        // advance layer iteration
        ++state.navigation.navLayerIter;
        return states::LayerToLayer{};  // nothing to do for this layer
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::SurfaceToSurface&,
                          const events::Status&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      // did we hit the layer?
      const Surface* surface = state.navigation.navSurfaceIter->representation;

      auto targetStatus =
          stepper.updateSurfaceStatus(state.stepping, *surface, true);

      if (targetStatus == Intersection3D::Status::onSurface) {
        ACTS_DEBUG(
            "Surface hit, current surface is now: " << surface->geometryId());
        state.navigation.currentSurface = surface;
        return std::nullopt;  // stay in state
      } else {
        ACTS_DEBUG("Surface not hit. Continue");
        state.navigation.currentSurface = nullptr;
        return std::nullopt;  // stay in state
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::SurfaceToSurface&,
                          const events::Target&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      // did we hit the surface?
      auto& ix = *state.navigation.navSurfaceIter;
      // auto& [init_ix, obj, iterSurface, alt] =
      // *state.navigation.navSurfaceIter;
      if (ix.representation == state.navigation.currentSurface) {
        // we're on the surface, just go on
        ++state.navigation.navSurfaceIter;
      }

      NavigationOptions<Surface> navOpts(state.stepping.navDir, true);

      while (state.navigation.navSurfaceIter !=
             state.navigation.navSurfaces.end()) {
        const Surface* surface =
            state.navigation.navSurfaceIter->representation;

        auto surfaceStatus =
            stepper.updateSurfaceStatus(state.stepping, *surface, true);
        if (surfaceStatus != Intersection3D::Status::reachable) {
          ACTS_DEBUG("During approach of surface "
                     << surface->geometryId()
                     << " intersection became unreachable => skipping surface "
                        "candidate");
          ++state.navigation.navSurfaceIter;
          continue;
        } else {
          // update straight line estimation
          ACTS_DEBUG("Proceeding towards surface: "
                     << surface->geometryId() << ", updated step size: "
                     << stepper.outputStepSize(state.stepping));
          return std::nullopt;  // stay in state
        }
      }

      if (state.navigation.navSurfaceIter ==
          state.navigation.navSurfaces.end()) {
        ACTS_DEBUG("No further surfaces in layer, back to layer to layer");
        ++state.navigation.navLayerIter;
        // state.navigation.navLayers.clear();
        // state.navigation.navLayerIter = stat.navigation.navLayers.end();
        return states::LayerToLayer{};
      }

      // we shouldn't get here
      return Terminated{};
    }

    template <typename propagator_state_t, typename stepper_t>
    void on_enter(const states::ToBoundarySurface&,
                  const FSMNavigator& navigator, propagator_state_t& state,
                  const stepper_t& stepper) {
      ACTS_DEBUG(volstr(*state.navigation.currentVolume)
                 << ": Targeting boundary surfaces");

      dispatch(events::ResolveBoundarySurfaces{}, navigator, state, stepper);
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::ToBoundarySurface&,
                          const events::ResolveBoundarySurfaces&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      ACTS_DEBUG(volstr(*state.navigation.currentVolume)
                 << ": Resolving boundary surfaces");

      // The navigation options
      NavigationOptions<Surface> navOpts(state.stepping.navDir, true);
      navOpts.pathLimit =
          state.stepping.stepSize.value(ConstrainedStep::aborter);
      navOpts.overstepLimit = stepper.overstepLimit(state.stepping);

      state.navigation.navBoundaries =
          state.navigation.currentVolume->compatibleBoundaries(
              state.geoContext, stepper.position(state.stepping),
              stepper.direction(state.stepping), navOpts,
              LoggerWrapper{logger()});

      // The number of boundary candidates
      if (logger().doPrint(Logging::VERBOSE)) {
        std::ostringstream os;
        os << state.navigation.navBoundaries.size();
        os << " boundary candidates found at path(s): ";
        for (auto& bc : state.navigation.navBoundaries) {
          os << bc.intersection.pathLength << "  ";
        }
        logger().log(Logging::VERBOSE, os.str());
      }
      // Set the begin iterator
      state.navigation.navBoundaryIter = state.navigation.navBoundaries.begin();
      if (not state.navigation.navBoundaries.empty()) {
        // Set to the first and return to the stepper
        stepper.updateStepSize(state.stepping,
                               *state.navigation.navBoundaryIter, true);
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": Navigation stepSize updated to "
                     << stepper.outputStepSize(state.stepping));

        return std::nullopt;
      } else {
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": No further bounding surfaces, is this bad?");
        return Terminated{};
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::ToBoundarySurface&,
                          const events::Status&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      ACTS_VERBOSE("ToBoundarySurface STATUS");
      // did we hit the boundary?
      const Surface* boundarySurface =
          state.navigation.navBoundaryIter->representation;

      auto targetStatus =
          stepper.updateSurfaceStatus(state.stepping, *boundarySurface, true);

      // @TODO: Check if we are in target volume

      if (targetStatus == Intersection3D::Status::onSurface) {
        ACTS_DEBUG("Boundary surface hit, current surface is now: "
                   << boundarySurface->geometryId());
        state.navigation.currentSurface = boundarySurface;
        return std::nullopt;  // stay in state
      } else {
        ACTS_DEBUG("Boundary not hit. Continue");
        // @TODO: Implement continuing if the first boundary is not hit
        state.navigation.currentSurface = nullptr;
        return std::nullopt;  // stay in state
      }
    }

    template <typename propagator_state_t, typename stepper_t>
    event_return on_event(const states::ToBoundarySurface&,
                          const events::Target&,
                          const FSMNavigator& /*navigator*/,
                          propagator_state_t& state, const stepper_t& stepper) {
      ACTS_VERBOSE("ToBoundarySurface TARGET");
      auto& [init_ix, boundary, surface, dir] =
          *state.navigation.navBoundaryIter;

      if (state.navigation.currentSurface == surface) {
        ACTS_DEBUG("Transition to volume connected by boundary surface")
        auto* nextVolume = boundary->attachedVolume(
            state.geoContext, stepper.position(state.stepping),
            stepper.direction(state.stepping), state.stepping.navDir);

        if (nextVolume == nullptr) {
          ACTS_VERBOSE(
              volstr(*state.navigation.currentVolume)
              << ": No more volume to progress to, stopping navigation.");
          state.navigation.currentVolume = nullptr;
          state.navigation.navigationBreak = true;
          stepper.releaseStepSize(state.stepping);
          return states::Finished{};  // done, the standard EndOfWorldReached
                                      // aborter should terminate
        }

        state.navigation.currentVolume = nextVolume;
        ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                     << ": Volume updated.");

        return states::VolumeToVolume{};
      }

      // not there yet, (re)set target

      // Loop over the boundary surface
      while (state.navigation.navBoundaryIter !=
             state.navigation.navBoundaries.end()) {
        // That is the current boundary surface
        auto boundarySurface = state.navigation.navBoundaryIter->representation;
        // Step towards the boundary surfrace
        auto boundaryStatus =
            stepper.updateSurfaceStatus(state.stepping, *boundarySurface, true);
        if (boundaryStatus == Intersection3D::Status::reachable) {
          ACTS_VERBOSE(volstr(*state.navigation.currentVolume)
                       << ": Boundary reachable, step size updated to "
                       << stepper.outputStepSize(state.stepping));
          return std::nullopt;
        } else {
          if (logger().doPrint(Logging::VERBOSE)) {
            std::ostringstream os;
            os << "Boundary ";
            os << std::distance(state.navigation.navBoundaryIter,
                                state.navigation.navBoundaries.end());
            os << " out of " << state.navigation.navBoundaries.size();
            os << " not reachable anymore, switching to next.";
            logger().log(Logging::VERBOSE, os.str());
          }

          // Increase the iterator to the next one
          ++state.navigation.navBoundaryIter;
          return std::nullopt;
        }
      }

      // if (state.navigation.navBoundaryIter ==
      //     state.navigation.navBoundaries.end()) {
      ACTS_ERROR("Failed to reach a boundary surface");
      return Terminated{};
      // }
    }

    // terminate enter handler: throw exception
    template <typename... Args>
    void on_enter(const Terminated&, Args&&...) {
      throw std::runtime_error("FSMNavigator FSM entered Terminated state");
    }

    // Logging methods

    template <typename S, typename E, typename S2>
    void on_process(const S&, const E&, const S2&) {
      ACTS_VERBOSE("FSM: transition: [" << S::name << "] + <" << E::name
                                        << "> = " << S2::name);
    }

    template <typename S, typename E>
    void on_process(const S&, const E&) {
      ACTS_VERBOSE("FSM: internal transition: [" << S::name << "] + <"
                                                 << E::name << ">");
    }

    template <typename E>
    void on_process(const E&) {
      ACTS_VERBOSE("FSM: process_event: <" << E::name << ">");
    }

    // catch all handlers

    template <typename S, typename E, typename... Args>
    event_return on_event(const S&, const E&, Args&&...) {
      return Terminated{};
    }

    template <typename S, typename... Args>
    void on_enter(const S&, Args&&...) {}

    template <typename S, typename... Args>
    void on_exit(const S&, Args&&...) {}

    const Logger& logger() const { return *m_logger; }

    const Logger* m_logger;
  };

  using state_type = State;

  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// The tolerance used to defined "reached"
    double tolerance = s_onSurfaceTolerance;

    /// Configuration for this Navigator
    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;
  };

  FSMNavigator(Config cfg, std::unique_ptr<const Logger> logger =
                               getDefaultLogger("FSMNavigator", Logging::INFO))
      : m_cfg(std::move(cfg)), m_logger(std::move(logger)){};

  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& state, const stepper_t& stepper) const {
    // ACTS_VERBOSE("Status call");
    state.navigation.m_logger = m_logger.get();
    state.navigation.dispatch(events::Status{}, *this, state, stepper);
  }

  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    // ACTS_VERBOSE("Target call");
    state.navigation.m_logger = m_logger.get();
    state.navigation.dispatch(events::Target{}, *this, state, stepper);
  }

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
