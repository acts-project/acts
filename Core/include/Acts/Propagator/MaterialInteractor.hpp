// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Material/MaterialProperties.hpp"
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
    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }
    // Do nothing if nothing is what is requested.
    if (not(multipleScattering or energyLoss or recordInteractions)) {
      return;
    }
    // We only have material interactions if there is potential material
    const auto* surface = state.navigation.currentSurface;
    if (not(surface and surface->surfaceMaterial())) {
      return;
    }

    // Determine how the update should be handled
    MaterialUpdateStage updateStage = MaterialUpdateStage::fullUpdate;
    if (surface == state.navigation.startSurface) {
      // We are at the start surface
      debugLog(state, [&] {
        return std::string("Update on start surface: post-update mode.");
      });
      updateStage = MaterialUpdateStage::postUpdate;
      // Or is it the target surface ?
    } else if (surface == state.navigation.targetSurface) {
      debugLog(state, [&] {
        return std::string("Update on target surface: pre-update mode");
      });
      updateStage = MaterialUpdateStage::preUpdate;
    } else {
      debugLog(state, [&] {
        return std::string("Update while pass through: full mode.");
      });
    }

    // Prepare relevant input particle properties
    const auto pos = stepper.position(state.stepping);
    const auto time = stepper.time(state.stepping);
    const auto dir = stepper.direction(state.stepping);
    const auto momentum = stepper.momentum(state.stepping);
    const auto q = stepper.charge(state.stepping);
    const auto qOverP = q / momentum;
    const auto mass = state.options.mass;
    const auto pdg = state.options.absPdgCode;
    const auto nav = state.stepping.navDir;
    const auto performCovarianceTransport = state.stepping.covTransport;

    // Determine the effective traversed material and its properties
    MaterialProperties slab =
        surface->surfaceMaterial()->materialProperties(pos, nav, updateStage);
    // Material exists but it's not real, i.e. vacuum; there is nothing to do
    if (not slab) {
      return;
    }

    // Correct the material properties for non-zero incidence
    const auto pathCorrection =
        surface->pathCorrection(state.geoContext, pos, dir);
    slab.scaleThickness(pathCorrection);

    // Compute contributions from interactions
    double varPhi = 0.0;
    double varTheta = 0.0;
    double varQOverP = 0.0;
    double Eloss = 0.0;
    if (multipleScattering and performCovarianceTransport) {
      // TODO use momentum before or after energy loss in backward mode?
      const auto theta0 =
          computeMultipleScatteringTheta0(slab, pdg, mass, qOverP, q);
      // sigmaTheta = theta0
      varTheta = theta0 * theta0;
      // sigmaPhi = theta0 / sin(theta)
      const auto sigmaPhi = theta0 * (dir.norm() / dir.z());
      varPhi = sigmaPhi * sigmaPhi;
    }
    // TODO just ionisation loss or full energy loss?
    if (energyLoss and performCovarianceTransport) {
      const auto sigmaQOverP =
          computeEnergyLossLandauSigmaQOverP(slab, pdg, mass, qOverP, q);
      varQOverP = sigmaQOverP * sigmaQOverP;
    }
    if (energyLoss) {
      Eloss = computeEnergyLossBethe(slab, pdg, mass, qOverP, q);
      debugLog(state, [=] {
        using namespace UnitLiterals;
        std::stringstream dstream;
        dstream << slab;
        dstream << " pdg=" << pdg;
        dstream << " mass=" << mass / 1_MeV << "MeV";
        dstream << " momentum=" << momentum / 1_GeV << "GeV";
        dstream << " energyloss=" << Eloss / 1_MeV << "MeV";
        return dstream.str();
      });
    }

    // To integrate process noise, we need to transport
    // the covariance to the current position in space
    // the 'true' indicates re-initializaiton of the further transport
    if (performCovarianceTransport) {
      stepper.covarianceTransport(state.stepping, true);
    }

    // update track parameters and covariance
    if (nav == NavigationDirection::forward) {
      // in forward propagation, energy decreases and variances increase
      const auto nextE = std::sqrt(mass * mass + momentum * momentum) - Eloss;
      // put particle at rest if energy loss is too large
      const auto nextP =
          (mass < nextE) ? std::sqrt(nextE * nextE - mass * mass) : 0;
      stepper.update(state.stepping, pos, dir, nextP, time);
      state.stepping.cov(ePHI, ePHI) += varPhi;
      state.stepping.cov(eTHETA, eTHETA) += varTheta;
      state.stepping.cov(eQOP, eQOP) += varQOverP;
    } else {
      // in backward propagation, energy increases and variances decreases
      const auto nextE = std::sqrt(mass * mass + momentum * momentum) + Eloss;
      const auto nextP = std::sqrt(nextE * nextE - mass * mass);
      stepper.update(state.stepping, pos, dir, nextP, time);
      // ensure variances stay positive even after noise removals
      const auto varPhiUp = state.stepping.cov(ePHI, ePHI) - varPhi;
      const auto varThetaUp = state.stepping.cov(eTHETA, eTHETA) - varTheta;
      const auto varQOverPUp = state.stepping.cov(eQOP, eQOP) - varQOverP;
      state.stepping.cov(ePHI, ePHI) = std::max(0.0, varPhiUp);
      state.stepping.cov(eTHETA, eTHETA) = std::max(0.0, varThetaUp);
      state.stepping.cov(eQOP, eQOP) = std::max(0.0, varQOverPUp);
    }

    result.materialInX0 += slab.thicknessInX0();
    result.materialInL0 += slab.thicknessInL0();
    // Record the interaction if requested
    if (recordInteractions) {
      MaterialInteraction mi;
      mi.position = pos;
      mi.time = time;
      mi.direction = dir;
      mi.deltaP = stepper.momentum(state.stepping) - momentum;
      mi.sigmaPhi2 = varPhi;
      mi.sigmaTheta2 = varTheta;
      mi.sigmaQoP2 = varQOverP;
      mi.surface = surface;
      mi.pathCorrection = pathCorrection;
      mi.materialProperties = slab;
      result.materialInteractions.push_back(std::move(mi));
    }
  }

  /// Material interaction has no pure observer.
  template <typename propagator_state_t>
  void operator()(propagator_state_t& /* unused */) const {}

 private:
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
