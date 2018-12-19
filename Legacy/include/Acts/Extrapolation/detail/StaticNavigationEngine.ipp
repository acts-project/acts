// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticNavigationEngine.ipp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Extrapolation/IMaterialEffectsEngine.hpp"
#include "Acts/Extrapolation/IPropagationEngine.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/BoundarySurfaceT.hpp"

template <class T>
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolveBoundaryT(
    Acts::ExtrapolationCell<T>& eCell,
    Acts::NavigationDirection   pDir) const
{
  EX_MSG_DEBUG(
      ++eCell.navigationStep,
      "navigation",
      "",
      "resolve boundary situation leaving '"
          << eCell.leadVolume->volumeName()
          << (int(pDir) > 0 ? "' along momentum." : "' opposite momentum."));
  // initialize the extrapolation code to progress
  ExtrapolationCode eCode = ExtrapolationCode::InProgress;
  // [1] ------------------------ fast boundary access : take straight line
  // estimates as navigation guide --------------
  auto boundaryIntersections
      = eCell.leadVolume->boundarySurfacesOrdered(*eCell.leadParameters, pDir);
  // some screen output
  EX_MSG_VERBOSE(
      eCell.navigationStep,
      "navigation",
      "",
      "found " << boundaryIntersections.size() << " boundary surfaces to try"
               << (eCell.onLastBoundary() ? " - starting from last boundary."
                                          : "."));
  // remember them for the slow acces
  std::map<const BoundarySurfaceT<TrackingVolume>*, bool> bSurfacesTried;

  for (auto& boundaryCandidate : boundaryIntersections) {
    // the surface of the
    const BoundarySurfaceT<TrackingVolume>* bSurfaceTV
        = boundaryCandidate.object;

    const Surface& bSurface = bSurfaceTV->surfaceRepresentation();
    // skip if it's the last boundary surface
    if (eCell.onLastBoundary() && &bSurface == eCell.lastBoundarySurface) {
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "navigation",
          bSurface.geoID().value(GeometryID::boundary_mask),
          "skipping this candidate boundary - identical to last boundary.");
      continue;
    }
    // check this boudnary, possible return codes are:
    // - SuccessPathLimit     : propagation to boundary caused PathLimit to be
    // fail @todo implement protection asainst far of tries
    // - SuccessMaterialLimit : boundary was reached and material update on
    // boundary reached limit
    // - InProgress           : boundary was reached and ready for continueing
    // the navigation
    // - UnSet                : boundary was not reached, try the next one
    // - FailureLoop          : next Volume was previous Volume
    eCode = handleBoundaryT<T>(eCell, *bSurfaceTV, pDir);
    CHECK_ECODE_SUCCESS(eCell, eCode);
    // Failure or Unset are not triggering a return, try more sophisticated
    // navigation
    if (!eCode.inProgress()) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "boundary surface not reached with " << eCode.toString()
                                                          << ", skipping.");
      // book keeping for the slow access not to try again the same stuff
      bSurfacesTried[bSurfaceTV] = false;
      // skip to the next surface if there's one
      continue;
    }
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   bSurface.geoID().value(GeometryID::boundary_mask),
                   "boundary surface handling yielded code "
                       << eCode.toString());
    // set that this was the last boundary surface
    eCell.lastBoundarySurface = &bSurface;
    // and return the code yielded by the handleBoundaryT
    return eCode;
  }
  // bail-out in case no boundary has been found for Fatras mode
  if (eCell.configurationMode(ExtrapolationMode::FATRAS)) {
    EX_MSG_VERBOSE(
        eCell.navigationStep, "navigation", "", "Fatras loop protection.");
    // get out of this, we will not simulate loopers
    return ExtrapolationCode::SuccessPathLimit;
  }

  // [2] ------------------------ slow boundary access : take all boundary
  // surfaces and simply try --------------
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "navigation",
                 "",
                 "fast boundary navigation did "
                 "not succeeed - trying slow "
                 "navigation now.");
  // ignore the ones you have tried already
  for (auto& bSurfaceTV : eCell.leadVolume->boundarySurfaces()) {
    // we tried this one already, no point to do it again
    if (bSurfacesTried.size()
        && bSurfacesTried.find(bSurfaceTV.get()) != bSurfacesTried.end()) {
      continue;
    }
    // skip if it's the last boundary surface
    const Surface& bSurface = bSurfaceTV->surfaceRepresentation();
    if (&bSurface == eCell.lastBoundarySurface) {
      continue;
    }
    // screen output
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   bSurface.geoID().value(GeometryID::boundary_mask),
                   "trying a boundary surface.");
    // there is now loop protection in the slow access, needs to be done by hand
    // check this boudnary, possible return codes are:
    // - SuccessPathLimit     : propagation to boundary caused PathLimit to be
    // fail @todo implement protection againg far of tries
    // - SuccessMaterialLimit : boundary was reached and material update on
    // boundary reached limit
    // - InProgress           : boundary was reached and ready for continueing
    // the navigation
    // - UnSet                : boundary was not reached, try the next one
    eCode = handleBoundaryT<T>(eCell, *bSurfaceTV, pDir);
    CHECK_ECODE_SUCCESS(eCell, eCode);
    // Failure or Unset are not triggering a return, try more sophisticated
    // navigation
    if (!eCode.inProgress()) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "boundary surface not reached with " << eCode.toString()
                                                          << ", skipping.");
      // skip to the next surface if there's one
      continue;
    }
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   "",
                   "boundary surface handling yielded code "
                       << eCode.toString());
    // set that this was the last boundary surface
    eCell.lastBoundarySurface = &bSurface;
    // and return the code yielded by the handleBoundaryT
    return eCode;
  }
  // [3] ------------------------ slowest boundary access : step-out-of-volume
  // approach -------------------------
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "navigation",
                 "",
                 "slow boundary "
                 "navigation did not "
                 "succeeed - trying "
                 "step-out-of-volume "
                 "approach now");
  for (auto& boundaryCandidate : boundaryIntersections) {
    // the surface of the
    const BoundarySurfaceT<TrackingVolume>* bSurfaceTV
        = boundaryCandidate.object;
    const Surface& bSurface = bSurfaceTV->surfaceRepresentation();
    // check this boudnary, possible return codes are:
    // - SuccessPathLimit     : propagation to boundary caused PathLimit to be
    // fail @todo implement protection againg far of tries
    // - SuccessMaterialLimit : boundary was reached and material update on
    // boundary reached limit
    // - InProgress           : boundary was reached and ready for continueing
    // the navigation
    // - UnSet                : boundary was not reached, try the next one
    // - FailureLoop          : next Volume was previous Volume
    eCode = handleBoundaryT<T>(eCell, *bSurfaceTV, pDir, true);
    CHECK_ECODE_SUCCESS(eCell, eCode);
    // Failure or Unset are not triggering a return, try more sophisticated
    // navigation
    if (!eCode.inProgress()) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "boundary surface not reached with " << eCode.toString()
                                                          << ", skipping.");
      // skip to the next surface if there's one
      continue;
    }
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   "",
                   "boundary surface handling yielded code "
                       << eCode.toString());
    // set that this was the last boundary surface
    eCell.lastBoundarySurface = &bSurfaceTV->surfaceRepresentation();
    // and return the code yielded by the handleBoundaryT
    return eCode;
  }
  // return it back
  EX_MSG_DEBUG(eCell.navigationStep,
               "navigation",
               "",
               "could not resolve the boundary situation. Exiting.");

  return ExtrapolationCode::FailureNavigation;
}

/** handle the failure - as configured */
template <class T>
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::handleBoundaryT(
    Acts::ExtrapolationCell<T>&                         eCell,
    const Acts::BoundarySurfaceT<Acts::TrackingVolume>& bSurfaceTV,
    Acts::NavigationDirection                           pDir,
    bool                                                stepout) const
{
  // get the bondary surface and compare with last one to prevent loops
  const Surface& bSurface = bSurfaceTV.surfaceRepresentation();
  // propagate the parameters to the boundary (force boundaryCheck to true in
  // case it is not a step-out trial), possible return codes :
  // - SuccessPathLimit : pathLimit reached during propagation
  // - InProgress       : boundary reached
  // - Recovered        : boundary not reached
  ExtrapolationCode eCode
      = m_cfg.propagationEngine->propagate(eCell,
                                           bSurface,
                                           pDir,
                                           {ExtrapolationMode::CollectBoundary},
                                           !stepout,
                                           eCell.destinationCurvilinear);
  CHECK_ECODE_SUCCESS(eCell, eCode);
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "navigation",
                 bSurface.geoID().value(GeometryID::boundary_mask),
                 "propagation with eCode " << eCode.toString());
  // check for progress
  if (eCode.inProgress()) {
    // check if the boundary solution is compatible with the radial direction of
    // the extrapolation
    if (!eCell.checkRadialCompatibility()) {
      // screen output for the radial compatibility check
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "radial compatbility check failed, radial direction is: "
                         << eCell.radialDirection);
      // if it's not jump back to the last valid lead parameters and return
      // ExtrapolationCode::Unset as a trigger
      eCell.leadParameters = eCell.lastLeadParameters;
      return ExtrapolationCode::Unset;
    }
    EX_MSG_VERBOSE(
        eCell.navigationStep,
        "navigation",
        bSurface.geoID().value(GeometryID::boundary_mask),
        "parameters on boundary surface created, moving to next volume.");
    // get the nextVolume - modify the position in case you have a step out
    // trial, take attachment otherwise
    const TrackingVolume* nextVolume = stepout
        ? m_cfg.trackingGeometry->lowestTrackingVolume(
              Vector3D(eCell.leadParameters->position()
                       + pDir * eCell.leadParameters->momentum().normalized()))
        : bSurfaceTV.attachedVolume(
              eCell.leadParameters->position(),
              eCell.leadParameters->momentum().normalized(),
              pDir);
    // check if we have no nextVolume : boundary rechaed @todo it's not really a
    // success
    if (!nextVolume) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "No next volume found of '"
                         << eCell.leadVolume->volumeName()
                         << "'. End of known world ?");

      // move the last step parameters to endParameters
      if (eCell.extrapolationSteps.size()) {
        // verbose screen message
        EX_MSG_VERBOSE(eCell.navigationStep,
                       "navigation",
                       bSurface.geoID().value(GeometryID::boundary_mask),
                       "parameters on surface turn into end parameters.");
        // the last step is transformed into endParameters
        auto& lStep         = eCell.extrapolationSteps.back();
        eCell.endParameters = std::move(lStep.parameters);
      }
      // return a successful reaching of the boundary
      return ExtrapolationCode::SuccessBoundaryReached;
    }
    // check if it is a boundary reached case
    // - geometrySignature change and configuration to stop then triggers a
    // Success
    bool stopAtThisBoundary
        = eCell.configurationMode(ExtrapolationMode::StopAtBoundary)
        && (nextVolume->geometrySignature()
            != eCell.leadVolume->geometrySignature());
    // fill the boundary into the cache if successfully hit boundary surface
    // loop protection - relaxed for the cases where you start from the boundary
    if (eCell.leadVolume == nextVolume) {
      // the start parameters where on the boundary already give a relaxed
      // return code
      if (&bSurface == eCell.lastBoundarySurface) {
        return ExtrapolationCode::Unset;
      }
      // give some screen output as of why this happens
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     " loop in volume '" << nextVolume->volumeName() << "'.")
      // return a loop failure, parameter deletion will be done by cache
      return ExtrapolationCode::FailureLoop;
    }
    // update the with the information of the layer material - will change the
    // leadParameters
    if (bSurface.associatedMaterial()) {
      // now handle the material, possible return codes:
      // - InProgress            : material update performed or not (depending
      // on material)
      // - SuccessMaterialLimit  : material limit reached & configured to stop
      // there
      eCode = m_cfg.materialEffectsEngine->handleMaterial(
          eCell, &bSurface, pDir, fullUpdate);
      CHECK_ECODE_SUCCESS(eCell, eCode);
    }
    // break if configured to break at volume boundary and signature change
    if (stopAtThisBoundary) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     bSurface.geoID().value(GeometryID::boundary_mask),
                     "geometry signature change from "
                         << eCell.leadVolume->geometrySignature()
                         << " to "
                         << nextVolume->geometrySignature());
      eCell.nextGeometrySignature = nextVolume->geometrySignature();
      // return the boundary reached : the navigation resolved already
      eCell.leadVolume = nextVolume;
      // move the last step parameters to endParameters
      if (eCell.extrapolationSteps.size()) {
        EX_MSG_VERBOSE(eCell.navigationStep,
                       "navigation",
                       bSurface.geoID().value(GeometryID::boundary_mask),
                       "parameters on surface turn into end parameters.");
        // the last step is transformed into endParameters
        auto& lStep         = eCell.extrapolationSteps.back();
        eCell.endParameters = std::move(lStep.parameters);
      }
      // give the return
      return ExtrapolationCode::SuccessBoundaryReached;
    }
    // remember the last boundary surface for loop protection
    eCell.lastBoundarySurface    = &bSurface;
    eCell.lastBoundaryParameters = eCell.leadParameters;
    // set next volume and reset lead layer
    eCell.leadVolume = nextVolume;
    eCell.leadLayer  = 0;
    // we have bParameters -> break the loop over boundaryIntersections
    return ExtrapolationCode::InProgress;
  }

  // you need to keep on trying
  return ExtrapolationCode::Unset;
}

/** handle the failure - as configured */
template <class T>
Acts::ExtrapolationCode
Acts::StaticNavigationEngine::resolvePositionT(
    Acts::ExtrapolationCell<T>& eCell,
    Acts::NavigationDirection   pDir,
    bool /*noLoop*/) const
{
  EX_MSG_DEBUG(
      ++eCell.navigationStep,
      "navigation",
      "",
      "resolve position (" << eCell.leadParameters->position().x() << ", "
                           << eCell.leadParameters->position().y()
                           << ", "
                           << eCell.leadParameters->position().z()
                           << ", "
                           << (int(pDir) > 0 ? ") along momentum."
                                             : ") opposite momentum."));

  return ExtrapolationCode::FailureNavigation;
}
