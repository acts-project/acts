// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <sstream>

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Propagator/detail/VolumeMaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

/// @brief The Material interaction struct
/// It records the surface  and the passed material
/// This is only nessecary recorded when configured
struct MaterialInteraction {
  /// The particle position at the interaction.
  Vector3D position = Vector3D(0., 0., 0);
  /// The particle time at the interaction.
  double time = 0.0;
  /// The particle direction at the interaction.
  Vector3D direction = Vector3D(0., 0., 0);
  /// The momentum change due to the interaction.
  double deltaP = 0.0;
  /// Expected phi variance due to the interactions.
  double sigmaPhi2 = 0.0;
  /// Expected theta variance due to the interactions.
  double sigmaTheta2 = 0.0;
  /// Expected q/p variance due to the interactions.
  double sigmaQoP2 = 0.0;
  /// The surface where the interaction occured.
  const Surface* surface = nullptr;
  /// The volume where the interaction occured.
  const TrackingVolume* volume = nullptr;
  /// Update the volume step to implment the proper step size
  bool updatedVolumeStep = false;
  /// The path correction factor due to non-zero incidence on the surface.
  /// The path correction factor due to non-zero incidence on the surface.
  double pathCorrection = 1.;
  /// The effective, passed material properties including the path correction.
  MaterialProperties materialProperties;
};

/// Material interactor propagator action.
///
/// Apply material interactions at a surface and update the track state.
struct MaterialInteractor {
  /// Whether to consider multiple scattering.
  bool multipleScattering = true;
  /// Whether to consider energy loss.
  bool energyLoss = true;
  /// Whether to record all material interactions.
  bool recordInteractions = false;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state which is
  /// created for every propagation/extrapolation step
  struct Result {
    // The accumulated materialInX0
    double materialInX0 = 0.;
    /// The accumulated materialInL0
    double materialInL0 = 0.;
    /// This one is only filled when recordInteractions is switched on
    std::vector<MaterialInteraction> materialInteractions;
  };
  using result_type = Result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // In case of Volume material update the result of the previous step
    if (recordInteractions && !result.materialInteractions.empty() &&
        result.materialInteractions.back().volume != nullptr &&
        result.materialInteractions.back().updatedVolumeStep == false) {
      UpdateResult(state, stepper, result);
    }

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // Do nothing if nothing is what is requested.
    if (not(multipleScattering or energyLoss or recordInteractions)) {
      return;
    }
    // We only have material interactions if there is potential material
    const Surface* surface = state.navigation.currentSurface;
    const TrackingVolume* volume = state.navigation.currentVolume;

    if (not(surface and surface->surfaceMaterial()) and
        not(volume and volume->volumeMaterial())) {
      return;
    }

    if (surface and surface->surfaceMaterial()) {
      // Prepare relevant input particle properties
      detail::PointwiseMaterialInteraction d(surface, state, stepper);

      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (not d.evaluateMaterialProperties(state)) {
        return;
      }

      // Evaluate the material effects
      d.evaluatePointwiseMaterialInteraction(multipleScattering, energyLoss);

      if (energyLoss) {
        debugLog(state, [=] {
          using namespace UnitLiterals;
          std::stringstream dstream;
          dstream << d.slab;
          dstream << " pdg=" << d.pdg;
          dstream << " mass=" << d.mass / 1_MeV << "MeV";
          dstream << " momentum=" << d.momentum / 1_GeV << "GeV";
          dstream << " energyloss=" << d.Eloss / 1_MeV << "MeV";
          return dstream.str();
        });
      }

      // To integrate process noise, we need to transport
      // the covariance to the current position in space
      if (d.performCovarianceTransport) {
        stepper.covarianceTransport(state.stepping);
      }
      // Apply the material interactions
      d.updateState(state, stepper);
      // Record the result
      recordResult(d, result);
    } else if (recordInteractions && volume and volume->volumeMaterial()) {
      // Prepare relevant input particle properties
      detail::VolumeMaterialInteraction d(volume, state, stepper);
      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (not d.evaluateMaterialProperties(state)) {
        return;
      }
      // Record the result
      recordResult(d, result);
    }
  }

  /// Material interaction has no pure observer.
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /* unused */) const {}

 private:
  /// @brief This function records the material effect
  ///
  /// @param [in] d Data cache container
  /// @param [in, out] result Result storage
  void recordResult(const detail::PointwiseMaterialInteraction& d,
                    result_type& result) const {
    result.materialInX0 += d.slab.thicknessInX0();
    result.materialInL0 += d.slab.thicknessInL0();
    // Record the interaction if requested
    if (recordInteractions) {
      MaterialInteraction mi;
      mi.position = d.pos;
      mi.time = d.time;
      mi.direction = d.dir;
      mi.deltaP = d.nextP - d.momentum;
      mi.sigmaPhi2 = d.variancePhi;
      mi.sigmaTheta2 = d.varianceTheta;
      mi.sigmaQoP2 = d.varianceQoverP;
      mi.surface = d.surface;
      mi.volume = nullptr;
      mi.pathCorrection = d.pathCorrection;
      mi.materialProperties = d.slab;
      result.materialInteractions.push_back(std::move(mi));
    }
  }

  /// @brief This function records the material effect
  ///
  /// @param [in] d Data cache container
  /// @param [in, out] result Result storage
  void recordResult(const detail::VolumeMaterialInteraction& d,
                    result_type& result) const {
    // Record the interaction
    MaterialInteraction mi;
    mi.position = d.pos;
    mi.time = d.time;
    mi.direction = d.dir;
    mi.sigmaPhi2 = d.variancePhi;
    mi.sigmaTheta2 = d.varianceTheta;
    mi.sigmaQoP2 = d.varianceQoverP;
    mi.surface = nullptr;
    mi.volume = d.volume;
    mi.pathCorrection = d.pathCorrection;
    mi.materialProperties = d.slab;
    result.materialInteractions.push_back(std::move(mi));
  }

  /// @brief This function update the previous material step
  ///
  /// @param [in] d Data cache container
  /// @param [in, out] result Result storage
  template <typename propagator_state_t, typename stepper_t>
  void UpdateResult(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
    // Update the previous interaction
    Vector3D shift = stepper.position(state.stepping) -
                     result.materialInteractions.back().position;
    double momentum = stepper.direction(state.stepping).norm();
    result.materialInteractions.back().deltaP =
        momentum - result.materialInteractions.back().direction.norm();
    result.materialInteractions.back().materialProperties.scaleThickness(
        shift.norm());
    result.materialInteractions.back().updatedVolumeStep = true;
    result.materialInX0 +=
        result.materialInteractions.back().materialProperties.thicknessInX0();
    result.materialInL0 +=
        result.materialInteractions.back().materialProperties.thicknessInL0();
  }

  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug ==
  /// true case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param state the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a streamable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(state.options.debugPfxWidth);
      dstream << "material interaction"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};

/// Using some short hands for Recorded Material
using RecordedMaterial = MaterialInteractor::Result;

/// And recorded material track
/// - this is start:  position, start momentum
///   and the Recorded material
using RecordedMaterialTrack =
    std::pair<std::pair<Acts::Vector3D, Acts::Vector3D>, RecordedMaterial>;

}  // namespace Acts
