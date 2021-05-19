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
  };

 public:
  struct State : FiniteStateMachine<State, states::Init, states::VolumeToVolume,
                                    states::LayerToLayer> {
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
    // const TrackingVolume* targetVolume = nullptr;
    /// Navigation state: the target layer
    // const Layer* targetLayer = nullptr;
    /// Navigation state: the target surface
    const Surface* targetSurface = nullptr;

    std::vector<LayerIntersection> navLayers;
    std::vector<LayerIntersection>::iterator navLayerIter;

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
      NavigationOptions<Layer> navOpts(
          state.stepping.navDir, true, navigator.m_cfg.resolveSensitive,
          navigator.m_cfg.resolveMaterial, navigator.m_cfg.resolvePassive,
          startLayer, nullptr);
      navOpts.pathLimit =
          state.stepping.stepSize.value(ConstrainedStep::aborter);

      state.navigation.navLayers =
          state.navigation.currentVolume->compatibleLayers(
              state.geoContext, stepper.position(state.stepping),
              stepper.direction(state.stepping), navOpts);

      ACTS_DEBUG("Found " << state.navigation.navLayers.size()
                          << " layer candidates");

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

      // updateStep(state, navCorr,
      //            state.navigation.navLayerIter->intersection.pathLength);

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
        state.navigation.currentSurface = layerSurface;
        if (layerSurface != nullptr) {
          ACTS_DEBUG("Layer hit, current surface is now"
                     << layerSurface->geometryId());
        }
        return std::nullopt;  // stay in state
      } else {
        ACTS_DEBUG("Layer not hit. Continue");
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
        // hm?
      } else {
        NavigationOptions<Surface> navOpts(state.stepping.navDir, true);

        while (state.navigation.navLayerIter !=
               state.navigation.navLayers.end()) {
          const Surface* layerSurface =
              state.navigation.navLayerIter->representation;

          auto layerStatus =
              stepper.updateSurfaceStatus(state.stepping, *layerSurface, true);
          if (layerStatus != Intersection3D::Status::reachable) {
            ACTS_DEBUG("During approach of layer" << layer->geometryId()
                                                  << " intersection");
            ACTS_DEBUG("became unreachable => skipping layer candidate");
            ++state.navigation.navLayerIter;
            continue;
          } else {
            // update straight line estimation
            ACTS_DEBUG("Proceeding towards layer: "
                       << layerSurface->geometryId() << ", updated step size: "
                       << stepper.outputStepSize(state.stepping));
            // updateStep(state, navCorr, ix.intersection.pathLength);
            break;
          }
        }
      }
      return std::nullopt;  // stay in state
    }

    template <typename propagator_state_t, typename stepper_t>
    void on_enter(const states::LayerToLayer&,
                  const FSMNavigator& /*navigator*/,
                  propagator_state_t& /*state*/, const stepper_t& /*stepper*/) {
      // ACTS_VERBOSE("Entering LayerToLayer loop: resolving layers!");
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
    std::shared_ptr<const TrackingGeometry> trackingGeometry;

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

 private:
  Config m_cfg;
  std::unique_ptr<const Logger> m_logger;

  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts
