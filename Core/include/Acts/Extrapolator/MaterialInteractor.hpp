// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Material/SurfaceMaterial.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

/// @brief The Material interaction struct
/// It records the surface  and the passed material
/// This is only nessecary recorded when configured
struct MaterialInteraction
{
  /// The material surface
  const Surface* surface = nullptr;

  /// The passsed materials
  /// it is the material and the actual (corrected) path length
  std::pair<Material, double> passedMaterial;

  double sigmaPhi2   = 0.;  ///< applied material effect: sigma(phi)^2
  double sigmaTheta2 = 0.;  ///< applied material effect: sigma(theta)^2
  double deltaP      = 0.;  ///< applied material effect: dela(p)
  double sigmaQoP2   = 0.;  ///< applied material effect: sigma(qop)^2
};

/// The Material interactor struct
///
/// This is a plugin to the Propagator that
/// performs material interaction on the currentSurface
/// of the Propagagor state
struct MaterialInteractor
{

  // Configuration for this MaterialInteractor
  /// multiple scattering switch on/off
  bool multipleScattering = true;
  /// the scattering struct
  detail::HighlandScattering scattering;

  /// energy loss switch on/off
  bool energyLoss = true;
  /// the energy loss struct
  detail::IonisationLoss ionisationloss;

  /// record material in detail
  bool recordDetailed = false;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state state which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    std::vector<MaterialInteraction> materialInteractions;
  };

  using result_type = this_result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  ///
  /// @param state is the mutable propagator state object
  /// @param result is the mutable result state object
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state, result_type& result) const
  {

    // if we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    // if switched off, then return - alows run-time configuration
    if (!multipleScattering && !energyLoss) {
      return;
    }

    // a current surface has been assigned by the navigator
    if (state.navigation.currentSurface
        && state.navigation.currentSurface->associatedMaterial()) {

      // get the surface material and the corresponding material properties
      auto sMaterial   = state.navigation.currentSurface->associatedMaterial();
      auto mProperties = sMaterial->material(state.stepping.position());
      if (mProperties) {
        // pre - full - post update test, i.e.
        // check if you have a factor for pre/post/full update to do
        double prepofu = 1.;
        if (state.navigation.startSurface == state.navigation.currentSurface) {
          debugLog(state, [&] {
            return std::string("Update on start surface: post-update mode.");
          });
          prepofu
              = state.navigation.currentSurface->associatedMaterial()->factor(
                  state.stepping.navDir, postUpdate);
        } else if (state.navigation.targetSurface
                   == state.navigation.currentSurface) {
          debugLog(state, [&] {
            return std::string("Update on target surface: pre-update mode");
          });
          prepofu
              = state.navigation.currentSurface->associatedMaterial()->factor(
                  state.stepping.navDir, preUpdate);
        } else {
          debugLog(state, [&] {
            return std::string("Update while pass through: full mode.");
          });
        }

        // create the material interaction class, in case we record afterwards
        MaterialInteraction mInteraction;
        mInteraction.surface = state.navigation.currentSurface;

        // the pre/post factor has been applied
        // now check if there's still something to do
        if (prepofu == 0.) {
          debugLog(state, [&] {
            return std::string("Pre/Post factor set material to zero.");
          });
          return;
        }

        // to integrate process noise, we need to transport
        // the covariance to the current position in space
        if (state.stepping.covTransport) {
          state.stepping.covarianceTransport(true);
        }

        // get the material thickness - and correct it with incidence
        double thickness = mProperties->thickness();
        // get the path correction due to the incident angle
        double pCorrection = state.navigation.currentSurface->pathCorrection(
            state.stepping.position(), state.stepping.direction());
        // the corrected thickness
        double cThickness = thickness * pCorrection;

        // the momentum at current position
        const double p     = state.stepping.p;
        const double m     = state.options.mass;
        const double E     = std::sqrt(p * p + m * m);
        const double lbeta = p / E;

        // apply the multiple scattering
        // - only when you do covariance transport
        if (multipleScattering && state.stepping.covTransport) {
          // thickness in X0 from without path correction
          double tInX0 = mProperties->thicknessInX0();
          // retrieve the scattering contribution
          double sigmaScat = scattering(p, lbeta, tInX0 * pCorrection);
          double sinTheta  = std::sin(state.stepping.direction().theta());
          double sigmaDeltaPhiSq
              = sigmaScat * sigmaScat / (sinTheta * sinTheta);
          double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
          // record the material interaction
          mInteraction.sigmaPhi2   = sigmaDeltaPhiSq;
          mInteraction.sigmaTheta2 = sigmaDeltaThetaSq;
          // good in any case for positive direction
          if (state.stepping.navDir == forward) {
            // just add the multiple scattering component
            state.stepping.cov(ePHI, ePHI)
                += state.stepping.navDir * sigmaDeltaPhiSq;
            state.stepping.cov(eTHETA, eTHETA)
                += state.stepping.navDir * sigmaDeltaThetaSq;
          } else {
            // we check if the covariance stays positive
            double sEphi   = state.stepping.cov(ePHI, ePHI);
            double sEtheta = state.stepping.cov(eTHETA, eTHETA);
            if (sEphi > sigmaDeltaPhiSq && sEtheta > sigmaDeltaThetaSq) {
              // noise removal is not applied if covariance would fall below 0
              state.stepping.cov(ePHI, ePHI) -= sigmaDeltaPhiSq;
              state.stepping.cov(eTHETA, eTHETA) -= sigmaDeltaThetaSq;
            }
          }
        }
        // apply the energy loss
        if (energyLoss) {
          // get the material
          const Material& mat = mProperties->material();
          // calculate gamma
          const double lgamma = E / m;
          // energy loss and straggling - per unit length
          std::pair<double, double> eLoss
              = ionisationloss(m, lbeta, lgamma, mat, 1. * units::_mm);
          // apply the energy loss
          const double dEdl = state.stepping.navDir * eLoss.first;
          const double dE   = thickness * pCorrection * dEdl;
          // check for energy conservation, and only apply momentum change
          // when kinematically allowed
          if (E + dE > m) {
            // calcuate the new momentum
            const double newP = std::sqrt((E + dE) * (E + dE) - m * m);
            // record the deltaP
            mInteraction.deltaP = p - newP;
            // update the state/momentum
            state.stepping.p = std::copysign(newP, state.stepping.p);
          }
          // transfer this into energy loss straggling and appply to covariance:
          // do that even if you had not applied energy loss do to
          // the kineamtic limit to catch the cases of deltE < MOP/MPV
          if (state.stepping.covTransport) {
            // calculate the straggling
            const double sigmaQoverP
                = thickness * pCorrection * eLoss.second / (lbeta * p * p);
            // save the material interaction
            mInteraction.sigmaQoP2 = sigmaQoverP * sigmaQoverP;
            // good in any case for positive direction
            if (state.stepping.navDir == forward) {
              state.stepping.cov(eQOP, eQOP)
                  += state.stepping.navDir * sigmaQoverP * sigmaQoverP;
            } else {
              // check that covariance entry doesn't become neagive
              double sEqop = state.stepping.cov(eQOP, eQOP);
              if (sEqop > sigmaQoverP * sigmaQoverP) {
                state.stepping.cov(eQOP, eQOP)
                    += state.stepping.navDir * sigmaQoverP * sigmaQoverP;
              }
            }
          }
        }
        // record if configured to do so
        if (recordDetailed) {
          // retrieves the material again (not optimal),
          // though this is not time critical
          const Material matr = mProperties->material();
          mInteraction.passedMaterial
              = std::pair<Material, double>(matr, cThickness);
          // record the material
          result.materialInteractions.push_back(mInteraction);
        }
      }
    }
  }

  /// Pure observer interface
  /// This does not apply to the surface collector
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& /*state*/) const
  {
  }

private:
  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param state the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
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
}
