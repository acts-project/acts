// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ExtrapolationEngine.ipp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Detector/TrackingGeometry.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include <iomanip>
#include <iostream>

template <class T>
Acts::ExtrapolationCode
Acts::ExtrapolationEngine::extrapolateT(Acts::ExtrapolationCell<T>& eCell,
                                        const Acts::Surface*        sf,
                                        Acts::PropDirection         dir,
                                        const Acts::BoundaryCheck& bcheck) const
{
  EX_MSG_DEBUG(eCell.navigationStep,
               "extrapolate",
               "",
               "starting extrapolation sequence.");
  // initialize the navigation
  ExtrapolationCode eCode = initNavigation<T>(eCell, sf, dir);
  EX_MSG_VERBOSE(
      eCell.navigationStep,
      "extrapolate",
      "",
      "initialize navigation with return code : " << eCode.toString());
  // main loop over volumes
  while (eCell.leadVolume && eCode == ExtrapolationCode::InProgress) {
    // give output that you are in the master volume loop
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "extrapolate",
                   "loop",
                   "processing volume : " << eCell.leadVolume->volumeName());
    // get the appropriate IExtrapolationEngine
    GeometryType geoType
        = eCell.leadVolume->geometrySignature() > 2 ? Dense : Static;

    std::shared_ptr<const IExtrapolationEngine> iee
        = (m_cfg.extrapolationEngines.size()
           > eCell.leadVolume->geometrySignature())
        ? m_cfg.extrapolationEngines[geoType]
        : nullptr;
    eCode = iee ? iee->extrapolate(eCell, sf, bcheck)
                : ExtrapolationCode::FailureConfiguration;
    // give a message about what you have
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "extrapolate",
                   "",
                   "returned from volume with return code : "
                       << eCode.toString()
                       << " and geoType:"
                       << geoType);
  }
  EX_MSG_DEBUG(
      eCell.navigationStep,
      "extrapolate",
      "",
      "extrapolation finished with return code : " << eCode.toString());
  // return the code
  return eCode;
}

template <class T>
Acts::ExtrapolationCode
Acts::ExtrapolationEngine::initNavigation(Acts::ExtrapolationCell<T>& eCell,
                                          const Acts::Surface*        sf,
                                          Acts::PropDirection         dir) const
{
  // initialize the Navigation stream
  // ----------------------------------------------------------------------------------------
  //
  // this is the global initialization, it only associated direct objects
  // detailed navigation search needs to be done by the sub engines (since they
  // know best)
  EX_MSG_DEBUG(++eCell.navigationStep,
               "navigation",
               "",
               "initialize the navigation stream.");
  // initialization of the navigation requires that leadParameters to be the
  // startParameters
  eCell.leadParameters = &(eCell.startParameters);
  // now check the tracking geometry and retrieve it if not existing
  if (!m_cfg.trackingGeometry) {
    EX_MSG_WARNING(eCell.navigationStep,
                   "navigation",
                   "",
                   "could not retrieve geometry. Stopping.");
    // configuration error
    return ExtrapolationCode::FailureConfiguration;
  } else
    EX_MSG_VERBOSE(
        eCell.navigationStep, "navigation", "", "geometry ready to use.");
  // ---------- START initialization
  // -----------------------------------------------------------------------------------------
  // initialize the start parameters - try association first
  eCell.startLayer = eCell.startLayer
      ? eCell.startLayer
      : eCell.leadParameters->associatedSurface().associatedLayer();
  eCell.startVolume = eCell.startVolume
      ? eCell.startVolume
      : (eCell.startLayer ? eCell.startLayer->enclosingTrackingVolume()
                          : nullptr);
  // check if you are at the volume boundary
  //!< @TODO do not use hard coded number of atVolumeBoundary, check with ST
  if (!eCell.startVolume
      || m_cfg.trackingGeometry->atVolumeBoundary(
             eCell.startParameters.position(), eCell.startVolume, 0.001)) {
    ExtrapolationCode ecVol
        = m_cfg.navigationEngine->resolvePosition(eCell, dir, true);
    if (!ecVol.isSuccessOrRecovered() && !ecVol.inProgress()) return ecVol;
    // the volume is found and assigned as start volume
    eCell.startVolume = eCell.leadVolume;
  } else {
    // we have a start volume, set is as lead volume also
    eCell.leadVolume = eCell.startVolume;
  }
  // bail out of the start volume can not be resolved
  if (!eCell.startVolume) return ExtrapolationCode::FailureNavigation;
  // screen output
  EX_MSG_VERBOSE(
      eCell.navigationStep,
      "navigation",
      "",
      "start volume termined as : " << eCell.startVolume->volumeName());
  // check layer association
  eCell.startLayer = eCell.startLayer
      ? eCell.startLayer
      : eCell.startVolume->associatedLayer(eCell.leadParameters->position());
  if (eCell.startLayer)
    EX_MSG_VERBOSE(
        eCell.navigationStep,
        "navigation",
        "",
        "start layer termined with index : " << eCell.startLayer->geoID());
  // the radial direction flags outwards or inwards moving directions w.r.t to
  // (0.,0.,0,)
  eCell.setRadialDirection();
  // ---------- END initialization
  // -----------------------------------------------------------------------------------------
  if (sf) {
    // stop at the end surface if configured
    eCell.endSurface = sf;
    // re-evaluate the radial direction if the end surface is given
    // should not happen in FATRAS extrapolation mode, usually Fatras has no end
    // surface though
    if (!eCell.checkConfigurationMode(ExtrapolationMode::FATRAS))
      eCell.setRadialDirection();
    // trying association via the layer : associated layer of material layer
    eCell.endLayer = sf->associatedLayer();
    eCell.endVolume
        = eCell.endLayer ? eCell.endLayer->enclosingTrackingVolume() : nullptr;
    // check if you found layer and volume
    if (!eCell.endVolume) {
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "navigation",
          "",
          "end volume needs to be determinded by surface intersection.");
      // use a propagation to find the endVolume and endLayer
      // @TODO can be opmisied (straight line for high momentum - use directly )
      ExtrapolationCell<T> navCell(*eCell.leadParameters, dir);
      // screen output
      ExtrapolationCode eCode = m_cfg.propagationEngine->propagate(
          navCell,
          *eCell.endSurface,
          anyDirection,
          ExtrapolationMode::Destination,
          false,
          eCell.navigationCurvilinear);
      // check for sucess to the destination
      //@TODO check what this is
      CHECK_ECODE_SUCCESS_NODEST(navCell, eCode);
      // screen output
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "navigation",
          "",
          "found endVolume and andLayer through propagation - return code : "
              << eCode.toString());
      // take the lead parameters to find end volume and end layer
      eCell.endVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          navCell.endParameters->position());
      eCell.endLayer = m_cfg.trackingGeometry->associatedLayer(
          navCell.endParameters->position());
    }
    // check the final end volume configuraiton - screen output
    if (eCell.endVolume)
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "navigation",
          "",
          "end volume termined as : " << eCell.endVolume->volumeName());
    if (eCell.endLayer)
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "navigation",
          "",
          "end layer termined with index : " << eCell.endLayer->geoID());
  } else
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   "",
                   "no destination surface nor end volume provided, "
                   "extrapolation has to stop by other means.");
  // return the progress call
  return ExtrapolationCode::InProgress;
}
