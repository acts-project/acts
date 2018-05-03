// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MATERIALINTERACTOR_H
#define ACTS_MATERIALINTERACTOR_H

#include <cmath>
#include <sstream>
#include <utility>
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/MaterialInteraction.hpp"

namespace Acts {

/// The Material interaction struct
/// - it records the surface
/// and the passed material
struct MaterialInteraction
{
  /// The material surface
  const Surface* surface = nullptr;
  /// The passsed materials
  /// it is the material and the actual (corrected) path length
  std::pair<Material, double> passedMaterial;
};

/// The Material interactor struct
///
/// @todo: check for covariance consistency
/// @todo: check for positive momentum
///
struct MaterialInteractor
{

  /// Configuration for this MaterialInteractor
  /// - multiple scattering
  bool multipleScattering = true;
  /// Configuration for this MaterialInteractor
  /// - energy loss
  bool energyLoss = true;
  /// use the mean (or most probable) energy loss
  bool energyLossMean = true;
  /// record material in detail
  bool recordDetailed = false;
  /// debug output flag
  bool debug = false;

  /// Simple result struct to be returned
  /// It mainly acts as an interal state cache which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    std::vector<MaterialInteraction> materialInteractions;
  };

  typedef this_result result_type;

  /// Interaction with detector material
  /// for the ActionList of the Propagator
  /// It checks if the cache has a current surface,
  /// in which case the action is performed:
  /// - the covariance is transported to the position,
  /// multiple scattering and energy loss is applied
  /// asccording to the configuration
  ///
  /// @tparam propagator_cache_t is the type of Propagagor cache
  /// @tparam stepper_cache_t is the type of Stepper cache
  ///
  /// @param pCache is the mutable propagator cache object
  /// @param sCache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& pCache,
             stepper_cache_t&    sCache,
             result_type&        result) const
  {

    // if we are on target, everything should have been done
    if (pCache.targetReached) return;

    // a current surface has been assigned by the navigator
    if (pCache.currentSurface && pCache.currentSurface->associatedMaterial()) {
      // @todo adapt sign & particle type
      ParticleType pType = muon;
      // get the surface material and the corresponding material properties
      auto sMaterial   = pCache.currentSurface->associatedMaterial();
      auto mProperties = sMaterial->material(sCache.position());
      if (mProperties) {
        // pre - full - post update test, i.e.
        // check if you have a factor for pre/post/full update to do
        double prepofu = 1.;
        if (pCache.startSurface == pCache.currentSurface) {
          debugLog(pCache, [&] {
            return std::string("Update on start surface: post-update mode.");
          });
          prepofu = pCache.currentSurface->associatedMaterial()->factor(
              sCache.navDir, postUpdate);
        } else if (pCache.targetSurface == pCache.currentSurface) {
          debugLog(pCache, [&] {
            return std::string("Update on target surface: pre-update mode");
          });
          prepofu = pCache.currentSurface->associatedMaterial()->factor(
              sCache.navDir, preUpdate);
        } else
          debugLog(pCache, [&] {
            return std::string("Update while pass through: full mode.");
          });
        // now check if there's still something to do
        if (prepofu == 0.) {
          debugLog(pCache, [&] {
            return std::string("Pre/Post factor set material to zero.");
          });
          return;
        }

        // to integrate process noise, we need to transport
        // the covariance to the current position in space
        sCache.applyCovTransport(true);
        // get the material thickness
        double thickness = mProperties->thickness();
        // get the path correction due to the incident angle
        double pCorrection = pCache.currentSurface->pathCorrection(
            sCache.position(), sCache.direction());
        // the corrected thickness
        double cThickness = thickness * pCorrection;
        // the momentum at current position
        double p    = std::abs(1. / sCache.qop);
        double m    = particleMasses.mass.at(pType);
        double E    = std::sqrt(p * p + m * m);
        double beta = p / E;
        // apply the multiple scattering
        if (multipleScattering) {
          // thickness in X0 from without path correction
          double tInX0 = mProperties->thicknessInX0();
          // retrieve the scattering contribution
          double sigmaScat = sigmaMS(tInX0 * pCorrection, p, beta);
          double sinTheta  = std::sin(sCache.direction().theta());
          double sigmaDeltaPhiSq
              = sigmaScat * sigmaScat / (sinTheta * sinTheta);
          double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
          // add or remove @todo implement check for covariance matrix -> 0
          sCache.cov(ePHI, ePHI) += sCache.navDir * sigmaDeltaPhiSq;
          sCache.cov(eTHETA, eTHETA) += sCache.navDir * sigmaDeltaThetaSq;
        }
        // apply the energy loss
        if (energyLoss) {
          // get the material
          const Material mat = mProperties->material();
          // energy loss and straggling
          std::pair<double, double> eLoss = energyLossMean
              ? ionizationEnergyLossMean(p, mat, pType)
              : ionizationEnergyLossMpv(p, mat, pType);
          // apply the energy loss
          const double dEdl   = sCache.navDir * eLoss.first;
          const double dE     = thickness * pCorrection * dEdl;
          double       sigmaP = eLoss.second;
          sigmaP *= thickness * pCorrection;
          // calcuate the new momentum
          const double newP = sqrt((E + dE) * (E + dE) - m * m);
          // update the cache/momentum
          std::copysign(1. / newP, sCache.qop);
          // transfer this into energy loss straggling
          const double sigmaDeltaE = thickness * pCorrection * sigmaP;
          const double sigmaQoverP = sigmaDeltaE / std::pow(beta * p, 2);
          // update the covariance if needed
          sCache.cov(eQOP, eQOP) += sCache.navDir * sigmaQoverP * sigmaQoverP;
        }
        // record if configured to do so
        if (recordDetailed) {
          // retrieves the material again (not optimal),
          // though this is not time critical
          const Material matr = mProperties->material();
          // create the material interaction class
          MaterialInteraction mInteraction;
          mInteraction.surface = pCache.currentSurface;
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
  template <typename propagator_cache_t, typename stepper_cache_t>
  void
  operator()(propagator_cache_t& pCache, stepper_cache_t& sCache) const
  {
    (void)pCache;
    (void)sCache;
  }

private:
  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the cache.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_cache_t Type of the propagator cache
  ///
  /// @param pCache the propagator cache for the debug flag, prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_cache_t>
  void
  debugLog(propagator_cache_t&          pCache,
           std::function<std::string()> logAction) const
  {
    if (debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(pCache.debugPfxWidth);
      dstream << "material interaction"
              << " | ";
      dstream << std::setw(pCache.debugMsgWidth) << logAction() << '\n';
      pCache.debugString += dstream.str();
    }
  }
};
}

#endif
