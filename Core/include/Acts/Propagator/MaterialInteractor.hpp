// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/PropagatorState.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/Propagator/detail/VolumeMaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

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
  /// Whether to add or remove noise.
  NoiseUpdateMode noiseUpdateMode = NoiseUpdateMode::addNoise;

  /// Type alias for material interaction result
  using result_type = RecordedMaterial;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper of the propagation
  /// @tparam navigator_t Type of the navigator of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param navigator The navigator in use
  /// @param result is the mutable result state object
  /// @param logger a logger instance
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                   const navigator_t& navigator, result_type& result,
                   const Logger& logger) const {
    if (state.stage == PropagatorStage::postPropagation) {
      return Result<void>::success();
    }

    // Do nothing if nothing is what is requested.
    if (!(multipleScattering || energyLoss || recordInteractions)) {
      return Result<void>::success();
    }

    // Handle surface material

    // Note that start and target surface conditions are handled in the
    // interaction code
    const Surface* surface = navigator.currentSurface(state.navigation);

    // We only have material interactions if there is potential material
    if (surface != nullptr && surface->surfaceMaterial()) {
      ACTS_VERBOSE("MaterialInteractor | " << "Found material on surface "
                                           << surface->geometryId());

      const MaterialSlab slab = detail::evaluateMaterialSlab(
          state, stepper, navigator,
          detail::determineMaterialUpdateMode(state, navigator,
                                              MaterialUpdateMode::FullUpdate));

      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (!slab.isVacuum()) {
        // To integrate process noise, we need to transport
        // the covariance to the current position in space
        if (state.stepping.covTransport) {
          stepper.transportCovarianceToCurvilinear(state.stepping);
        }

        const double initialMomentum = stepper.absoluteMomentum(state.stepping);

        // Apply the material interactions
        const detail::PointwiseMaterialEffects effects =
            detail::performMaterialInteraction(state, stepper, slab,
                                               noiseUpdateMode,
                                               multipleScattering, energyLoss);

        if (energyLoss) {
          using namespace UnitLiterals;

          const ParticleHypothesis& particleHypothesis =
              stepper.particleHypothesis(state.stepping);
          const PdgParticle absPdg = particleHypothesis.absolutePdg();
          const double mass = particleHypothesis.mass();
          const double momentum = stepper.absoluteMomentum(state.stepping);

          ACTS_VERBOSE("MaterialInteractor | "
                       << slab << " absPdg=" << absPdg
                       << " mass=" << mass / 1_MeV << "MeV"
                       << " momentum=" << momentum / 1_GeV << "GeV"
                       << " energyloss=" << effects.eLoss / 1_MeV << "MeV");
        }

        // Record the result
        recordResult(state, stepper, navigator, slab, initialMomentum, effects,
                     result);
      }
    }

    // Handle volume material

    // In case of Volume material update the result of the previous step
    if (!result.materialInteractions.empty() &&
        !result.materialInteractions.back().volume.empty() &&
        result.materialInteractions.back().updatedVolumeStep == false) {
      updateResult(state, stepper, result);
    }

    const TrackingVolume* volume = navigator.currentVolume(state.navigation);

    // We only have material interactions if there is potential material
    if (volume && volume->volumeMaterial()) {
      ACTS_VERBOSE("MaterialInteractor | " << "Found material in volume "
                                           << volume->geometryId());

      // Prepare relevant input particle properties
      detail::VolumeMaterialInteraction interaction(volume, state, stepper);
      // Determine the effective traversed material and its properties
      // Material exists but it's not real, i.e. vacuum; there is nothing to do
      if (interaction.evaluateMaterialSlab(state, navigator)) {
        // Record the result
        recordResult(interaction, result);
      }
    }

    return Result<void>::success();
  }

 private:
  /// @brief This function records the material effect
  ///
  /// @tparam propagator_state_t is the type of Propagator state
  /// @tparam stepper_t Type of the stepper of the propagation
  /// @tparam navigator_t Type of the navigator of the propagation
  ///
  /// @param [in] state The propagator state
  /// @param [in] stepper The stepper instance
  /// @param [in] navigator The navigator instance
  /// @param [in] slab The material slab
  /// @param [in] initialMomentum Initial momentum before the interaction
  /// @param [in] effects The material effects
  /// @param [in, out] result Result storage
  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void recordResult(const propagator_state_t& state, const stepper_t& stepper,
                    const navigator_t& navigator, const MaterialSlab& slab,
                    double initialMomentum,
                    const detail::PointwiseMaterialEffects& effects,
                    result_type& result) const {
    result.materialInX0 += slab.thicknessInX0();
    result.materialInL0 += slab.thicknessInL0();

    // Record the interaction if requested
    if (!recordInteractions) {
      return;
    }

    const Vector3 position = stepper.position(state.stepping);
    const double time = stepper.time(state.stepping);
    const Vector3 direction = stepper.direction(state.stepping);
    const double finalMomentum = stepper.absoluteMomentum(state.stepping);

    const Surface* surface = navigator.currentSurface(state.navigation);
    const double pathCorrection =
        surface->pathCorrection(state.options.geoContext, position, direction);

    MaterialInteraction mi;
    mi.position = position;
    mi.time = time;
    mi.direction = direction;
    mi.deltaP = initialMomentum - finalMomentum;
    mi.sigmaPhi2 = effects.variancePhi;
    mi.sigmaTheta2 = effects.varianceTheta;
    mi.sigmaQoP2 = effects.varianceQoverP;
    mi.surface = surface;
    mi.volume = InteractionVolume();
    mi.pathCorrection = pathCorrection;
    mi.materialSlab = slab;
    result.materialInteractions.push_back(std::move(mi));
  }

  /// @brief This function records the material effect
  ///
  /// @param [in] interaction Interaction cache container
  /// @param [in, out] result Result storage
  void recordResult(const detail::VolumeMaterialInteraction& interaction,
                    result_type& result) const {
    // Record the interaction if requested
    if (!recordInteractions) {
      return;
    }

    MaterialInteraction mi;
    mi.position = interaction.pos;
    mi.time = interaction.time;
    mi.direction = interaction.dir;
    mi.surface = nullptr;
    mi.volume = interaction.volume;
    mi.pathCorrection = interaction.pathCorrection;
    mi.materialSlab = interaction.slab;
    result.materialInteractions.push_back(std::move(mi));
  }

  /// @brief This function update the previous material step
  ///
  /// @param [in,out] state The state object
  /// @param [in] stepper The stepper instance
  /// @param [in, out] result Result storage
  template <typename propagator_state_t, typename stepper_t>
  void updateResult(propagator_state_t& state, const stepper_t& stepper,
                    result_type& result) const {
    // Update the previous interaction if requested
    if (!recordInteractions) {
      return;
    }

    Vector3 shift = stepper.position(state.stepping) -
                    result.materialInteractions.back().position;
    double momentum = stepper.direction(state.stepping).norm();
    result.materialInteractions.back().deltaP =
        momentum - result.materialInteractions.back().direction.norm();
    result.materialInteractions.back().materialSlab.scaleThickness(
        shift.norm());
    result.materialInteractions.back().updatedVolumeStep = true;
    result.materialInX0 +=
        result.materialInteractions.back().materialSlab.thicknessInX0();
    result.materialInL0 +=
        result.materialInteractions.back().materialSlab.thicknessInL0();
  }
};

}  // namespace Acts
