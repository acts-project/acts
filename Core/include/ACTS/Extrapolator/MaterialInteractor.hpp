// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MATERIALINTERACTOR_H
#define ACTS_MATERIALINTERACTOR_H

#include <sstream>
#include "ACTS/EventData/ParticleDefinitions.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Material/MaterialProperties.hpp"
#include "ACTS/Material/SurfaceMaterial.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/MaterialInteraction.hpp"

namespace Acts {

struct MaterialInteraction
{
  /// The material surface
  const Surface* surface = nullptr;
  /// The passsed materials
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
  bool multiple_scattering = true;
  /// Configuration for this MaterialInteractor
  /// - energy loss
  bool energy_loss = true;
  /// use the mean (or most probable) energy loss
  bool energy_loss_mean = true;
  /// record material in detail
  bool record_detailed = false;
  /// The particle masses for lookup
  ParticleMasses particleMasses;

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
  /// @tparam cache_t is the type of Stepper cache
  ///
  /// @param cache is the mutable stepper cache object
  /// @param result is the mutable result cache object
  template <typename cache_t>
  void
  operator()(cache_t& cache, result_type& result) const
  {
    // a current surface has been assigned by the navigator
    if (cache.current_surface && cache.current_surface->associatedMaterial()) {
      // @todo adapt sign & particle type
      int          sign  = 1;
      ParticleType pType = muon;
      // get the surface material and the corresponding material properties
      auto sMaterial   = cache.current_surface->associatedMaterial();
      auto mProperties = sMaterial->material(cache.position());
      if (mProperties) {
        // if we process noise, we need to transport
        // the covariance to the current place
        cache.apply_cov_transport(true);
        // get the material thickness
        double thickness = mProperties->thickness();
        // get the path correction due to the incident angle
        double pCorrection = cache.current_surface->pathCorrection(
            cache.position(), cache.direction());
        // the corrected thickness
        double cThickness = thickness * pCorrection;
        // the momentum at current position
        double p    = std::abs(1. / cache.qop);
        double m    = particleMasses.mass.at(pType);
        double E    = sqrt(p * p + m * m);
        double beta = p / E;
        // apply the multiple scattering
        if (multiple_scattering) {
          // thickness in X0 from without path correction
          double tInX0 = mProperties->thicknessInX0();
          // retrieve the scattering contribution
          double sigmaScat = sigmaMS(tInX0 * pCorrection, p, beta);
          double sinTheta  = sin(cache.direction().theta());
          double sigmaDeltaPhiSq
              = sigmaScat * sigmaScat / (sinTheta * sinTheta);
          double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
          // add or remove @todo implement check for covariance matrix -> 0
          cache.cov(ePHI, ePHI) += sign * sigmaDeltaPhiSq;
          cache.cov(eTHETA, eTHETA) += sign * sigmaDeltaThetaSq;
        }
        // apply the energy loss
        if (energy_loss) {
          // get the material
          const Material mat = mProperties->material();
          // energy loss and straggling
          std::pair<double, double> eLoss = energy_loss_mean
              ? ionizationEnergyLoss_mean(p, mat, pType, particleMasses)
              : ionizationEnergyLoss_mop(p, mat, pType, particleMasses);
          // apply the energy loss
          const double dEdl   = sign * eLoss.first;
          const double dE     = thickness * pCorrection * dEdl;
          double       sigmaP = eLoss.second;
          sigmaP *= thickness * pCorrection;
          // calcuate the new momentum
          const double newP = sqrt((E + dE) * (E + dE) - m * m);
          // update the cache/momentum
          double charge = cache.qop > 0. ? 1. : -1.;
          cache.qop     = charge / newP;
          // transfer this into energy loss straggling
          const double sigmaDeltaE = thickness * pCorrection * sigmaP;
          const double sigmaQoverP = sigmaDeltaE / std::pow(beta * p, 2);
          // update the covariance if needed
          cache.cov(eQOP, eQOP) += sign * sigmaQoverP * sigmaQoverP;
        }
        // record if configured to do so
        if (record_detailed) {
          // retrieves the material again (not optimal),
          // though this is not time critical
          const Material matr = mProperties->material();
          // create the material interaction class
          MaterialInteraction mInteraction;
          mInteraction.surface = cache.current_surface;
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
  template <typename cache_t>
  void
  operator()(cache_t& cache) const
  {
    (void)cache;
  }
};
}

#endif
