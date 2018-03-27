// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// StaticEngine.ipp, Acts project
///////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/Extrapolation/IMaterialEffectsEngine.hpp"
#include "Acts/Extrapolation/INavigationEngine.hpp"
#include "Acts/Extrapolation/IPropagationEngine.hpp"
#include "Acts/Layers/Layer.hpp"
#include "Acts/Surfaces/Surface.hpp"

template <class T>
Acts::ExtrapolationCode
Acts::StaticEngine::extrapolateT(Acts::ExtrapolationCell<T>& eCell,
                                 const Acts::Surface*        sf,
                                 Acts::NavigationDirection   pDir,
                                 const Acts::BoundaryCheck&  bcheck) const
{
  Acts::ExtrapolationCode eCode = Acts::ExtrapolationCode::InProgress;
  // ---- [0] check the direct propagation exit
  //
  //  obviously need a surface to exercise the fallback & need to be configured
  //  to do so
  if (sf && eCell.configurationMode(ExtrapolationMode::Destination)) {
    EX_MSG_DEBUG(
        ++eCell.navigationStep,
        "extrapolate",
        "",
        "direct extapolation in volume : " << eCell.leadVolume->volumeName());
    // propagate to the surface, possible return codes are : SuccessPathLimit,
    // SucessDestination, FailureDestination
    eCode = m_cfg.propagationEngine->propagate(eCell,
                                               *sf,
                                               pDir,
                                               {ExtrapolationMode::Destination},
                                               bcheck,
                                               eCell.destinationCurvilinear);
    // eCode can be directly returned
    return eCode;
  }

  EX_MSG_DEBUG(++eCell.navigationStep,
               "extrapolate",
               "",
               "extrapolation in static environment in volume : "
                   << eCell.leadVolume->volumeName());
  // evoke or finish the navigation initialization, possible return codes are:
  // - InProgress        : everything is fine, extrapolation in static volume is
  //                       in progress
  // - FailureNavigation : navigation setup could not be resolved, but reovery
  //                       was not configured
  // - Recovered         : navigation setup could not be resolved, recovery by
  //                       fallback to directly kicked in (and worked)
  eCode = initNavigationT<T>(eCell, sf, pDir, bcheck);
  CHECK_ECODE_CONTINUE(eCell, eCode);

  // ----- [1] handle the (leadLayer == endLayer ) case
  // this case does not need a layer to layer loop
  // - we directly resolve the layer and return
  if (sf && eCell.leadLayer == eCell.endLayer && eCell.initialVolume()) {
    // screen output for startLayer == endLayer
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   eCell.leadLayer->geoID().value(GeometryID::layer_mask),
                   "start and destination layer are identical.");
    // set the leadLayerSurface to the parameters surface for a start point
    eCell.leadLayerSurface = &(eCell.leadParameters->referenceSurface());
    // resolve the layer, it is the final extrapolation
    // - InProgress           : layer resolving went without problem
    // - SuccessPathLimit     : path limit reached & configured to stop
    //                          (rather unlikely within a layer)
    // - SuccessMaterialLimit : material limit successfully reached
    eCode = handleLayerT<T>(eCell, sf, pDir, bcheck);
    // Success triggers a return
    CHECK_ECODE_SUCCESS_NODEST(eCell, eCode);
    // extrapolation to destination was not successful
    // - handle the return as configured (e.g. fallback)
    return handleReturnT<T>(eCode, eCell, sf, pDir, bcheck);
  }

  // ----- [2] do the (layer-to-layer) loop
  //
  // the volume returns the layers ordered by distance
  // - start and end layer will be part of the loop
  // - the content of layers is given by the configuration mode
  // -- the surface on approach is already resolved
  // shortcuts for the collection bools
  bool collectSensitive
      = eCell.configurationMode(ExtrapolationMode::CollectSensitive);
  bool collectMaterial
      = eCell.configurationMode(ExtrapolationMode::CollectMaterial);
  bool collectPassive
      = eCell.configurationMode(ExtrapolationMode::CollectPassive);

  // check for the final volume
  const Layer* fLayer = eCell.finalVolumeReached() ? eCell.endLayer : nullptr;

  // get the layer intersections
  // this provides already the surfaceOnApproach when needed
  auto layerIntersections
      = eCell.leadVolume->layerCandidatesOrdered(eCell.leadLayer,
                                                 fLayer,
                                                 *eCell.leadParameters,
                                                 pDir,
                                                 true,
                                                 collectSensitive,
                                                 collectMaterial,
                                                 collectPassive);
  // some screen output
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "layer",
                 "loop",
                 "found " << layerIntersections.size()
                          << " layers for the layer-to-layer loop.");
  // the actual
  // layer-to-layer loop starts here
  for (auto& layerCandidate : layerIntersections) {
    // assign the leadLayer
    eCell.leadLayer = layerCandidate.object;
    // screen output for layer-to-layer loop
    EX_MSG_VERBOSE(
        eCell.navigationStep,
        "layer",
        "loop",
        "processing layer with index : "
            << eCell.leadLayer->geoID().value(GeometryID::layer_mask));
    // we set the lead layer surface :
    // - to the start surface for the start layer
    // - to the provided approach surface
    eCell.leadLayerSurface = (eCell.leadLayer == eCell.startLayer)
        ? &(eCell.leadParameters->referenceSurface())
        : layerCandidate.representation;
    // handle the layer, possible returns are :
    // - InProgress               : fine, whatever happened on the lead layer
    // - SuccessWithPathLimit     : propagation towards layer exceeded limit
    // - SuccessWithMaterialLimit : material interaction killed track
    // - FailureDestination       : destination was not hit appropriately
    eCode = handleLayerT<T>(eCell, sf, pDir, bcheck);
    // some more screen output
    EX_MSG_VERBOSE(
        eCell.navigationStep,
        "layer",
        layerCandidate.object->geoID().value(GeometryID::layer_mask),
        "handleLayerT returned extrapolation code : " << eCode.toString());
    // Possibilities are:
    // - SuccessX  -> return (via handleReturnT)
    // - FailureX  -> return (via handleReturnT that might evoke a fallback)
    // - InProgess -> continue layer-to-layer loop
    if (!eCode.inProgress())
      return handleReturnT<T>(eCode, eCell, sf, pDir, bcheck);
  }

  // ----- [3] now resolve the boundary situation, call includes information
  // wheather one is alreay at a boundary
  //
  // the navigaiton engine ca trigger different return codes
  // - InProgress                   : fine, boundary surface has been found
  // - SuccessWithPathLimit         : propagation towards boundary surface
  // exceeded path limit
  // - FailureLoop/Navigation       : problem in boundary resolving
  eCode = m_cfg.navigationEngine->resolveBoundary(eCell, pDir);
  // SuccessX and FailureX trigger a return
  CHECK_ECODE_SUCCESS_NODEST(eCell, eCode);
  // handle the return of the boudnary resolving
  return handleReturnT<T>(eCode, eCell, sf, pDir, bcheck);
}

template <class T>
Acts::ExtrapolationCode
Acts::StaticEngine::initNavigationT(Acts::ExtrapolationCell<T>& eCell,
                                    const Acts::Surface*        sf,
                                    Acts::NavigationDirection   pDir,
                                    const Acts::BoundaryCheck&  bcheck) const
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
               "complete for static environment.");
  // [A] the initial volume
  if (eCell.startVolume == eCell.leadVolume && eCell.startLayer) {
    // - found the initial start layer through association
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   "",
                   "this is the initial volume, everything set up already.");
    // assigning it to the leadLayer
    eCell.leadLayer = eCell.startLayer;
    // return progress
    return Acts::ExtrapolationCode::InProgress;
  }
  // [B] any volume if we don't have a leadLayer
  if (!eCell.leadLayer) {
    // - finding it through global search, never a boundary layer ... convention
    // says that you update by exit
    eCell.leadLayer
        = eCell.leadVolume->associatedLayer(eCell.leadParameters->position());
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "navigation",
                   "",
                   "no start layer found yet, looking for it ..."
                       << OH_CHECKFOUND(eCell.leadLayer));
  }
  // [C] the final volume - everything's fine
  if (eCell.leadVolume == eCell.endVolume && sf) {
    if (eCell.endLayer) {
      // the end layer had been found already by association
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "navigation",
                     "",
                     "this is the final volume, everything set up already.");
      return ExtrapolationCode::InProgress;
    } else {
      // make a straight line intersection
      Acts::Intersection sfI = sf->intersectionEstimate(
          eCell.leadParameters->position(),
          pDir * eCell.leadParameters->momentum().unit(),
          true);
      // use this to find endVolume and endLayer
      eCell.endLayer = eCell.leadVolume->associatedLayer(sfI.position);
      // if you have a surface you need to require an end layer for the
      // validation, otherwise you need to do a fallbac
      return eCell.endLayer
          ? Acts::ExtrapolationCode::InProgress
          : handleReturnT<T>(
                ExtrapolationCode::FailureNavigation, eCell, sf, pDir, bcheck);
    }
  }
  // return that you're in progress
  return Acts::ExtrapolationCode::InProgress;
}

template <class T>
Acts::ExtrapolationCode
Acts::StaticEngine::handleLayerT(ExtrapolationCell<T>& eCell,
                                 const Surface*        sf,
                                 NavigationDirection   pDir,
                                 const BoundaryCheck&  bcheck) const
{

  // initialize the extrapolaiton code
  ExtrapolationCode eCode = ExtrapolationCode::InProgress;

  /// prepare the layer output number
  auto layerValue = eCell.leadLayer->geoID().value(GeometryID::layer_mask);
  // screen output
  EX_MSG_DEBUG(++eCell.navigationStep,
               "layer",
               layerValue,
               "in volume " << eCell.leadVolume->volumeName());

  // start layer, i.e. no propagation needed
  bool isStartLayer
      = (eCell.leadLayerSurface == &eCell.leadParameters->referenceSurface());
  // final layer, end propagation needed
  bool isFinalLayer = (eCell.leadLayer == eCell.endLayer);

  // shortcuts for the collection bools
  bool collectSensitive
      = eCell.configurationMode(ExtrapolationMode::CollectSensitive);
  bool collectMaterial
      = eCell.configurationMode(ExtrapolationMode::CollectMaterial);
  bool collectPassive
      = eCell.configurationMode(ExtrapolationMode::CollectPassive);

  // ----- [0] the start situation on the layer needs to be resolved
  //
  // [A] the layer is not the start layer,
  //     - it needs propagation towards it on approach
  //     - the surface on approach has already been sorted out in the loop
  if (!isStartLayer) {
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   layerValue,
                   "not the start layer, propagate to it.");
    // propagate to the representing surface of this layer
    // - InProgress       : propagation to approaching surface worked
    //                      check material update
    // - SuccessPathLimit : propagation to approaching surface reached path
    // limit
    //                      limit
    // - Recovered        : layer was not hit, so can be ignored in the layer to
    //                      layer loop
    eCode = m_cfg.propagationEngine->propagate(
        eCell,
        *eCell.leadLayerSurface,
        pDir,
        {ExtrapolationMode::CollectPassive},
        true,
        eCell.navigationCurvilinear);
    CHECK_ECODE_SUCCESS_NODEST(eCell, eCode);
    // the extrapolation to the initial layer did not succeed
    // -> skip this layer and continue the layer-to-layer loop
    if (eCode == ExtrapolationCode::Recovered) {
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "layer",
                     layerValue,
                     "has not been hit, return to layer-to-layer loop.");
      return ExtrapolationCode::InProgress;
    }
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   layerValue,
                   "successfuly hit on approach");
    // the correct material layer needs to be assigned - in case of the approach
    // surface not being hit, his can be the layer surface
    if (eCell.leadLayerSurface->associatedMaterial()) {
      // now handle the material (full update when passing approach surface),
      // return codes are:
      // - SuccessMaterialLimit : material limit reached, return back
      // - InProgress           : material update done or not (depending on the
      // material description)
      eCode = m_cfg.materialEffectsEngine->handleMaterial(
          eCell, eCell.leadLayerSurface, pDir, fullUpdate);
      CHECK_ECODE_CONTINUE(eCell, eCode);
    }
  } else if (isStartLayer) {
    // [B] the layer is the start layer
    // - the parameters are already on the layer
    // - let's check if a post update on the start surface has to be done
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   layerValue,
                   "is start layer, no propagation to be done.");
    // the start surface could have a material layer attached
    const Surface& surface = eCell.leadParameters->referenceSurface();
    if (surface.associatedMaterial()) {
      // now handle the material (post update on start layer), return codes are:
      // - SuccessMaterialLimit : material limit reached, return back
      // - InProgress           : material update done or not
      //                          (depending on the material description)
      eCode = m_cfg.materialEffectsEngine->handleMaterial(
          eCell, &surface, pDir, postUpdate);
      CHECK_ECODE_CONTINUE(eCell, eCode);
    }
  }

  // ----- [1] the sub structure of the layer needs to be resolved:
  // this will give you the compatible surfaces of the layer
  // - provided start and destination surface are excluded
  //   the start surface already the one of the leadParameters
  // - surfaces without material are only provided if they are active and
  //   CollectSensitive is configured
  // - surfaces with material are provided in order to make the necessary
  //   material update
  //
  // the ordered surfaces are returned
  std::vector<Acts::SurfaceIntersection> cSurfaces
      = eCell.leadLayer->compatibleSurfaces(*eCell.leadParameters,
                                            pDir,
                                            bcheck,
                                            collectSensitive,
                                            collectMaterial,
                                            collectPassive,
                                            eCell.searchMode,
                                            eCell.leadLayerSurface,
                                            (isFinalLayer ? sf : nullptr));
  // how many test surfaces do we have
  size_t ncSurfaces = cSurfaces.size();

  // some screen output for the sub structure
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "layer",
                 layerValue,
                 "found " << ncSurfaces << " sub structure surfaces to test.");
  // check if you have to do something
  if (ncSurfaces) {

    // now loop over the surfaces - they are assumed to be sorted
    for (auto& csf : cSurfaces) {
      // get the surface
      auto         surface = csf.object;
      geo_id_value sensitiveID
          = surface->geoID().value(GeometryID::sensitive_mask);
      geo_id_value surfaceID = sensitiveID
          ? sensitiveID
          : surface->geoID().value(GeometryID::approach_mask);
      // the surface to try
      EX_MSG_VERBOSE(eCell.navigationStep,
                     "surface",
                     surfaceID,
                     "trying candidate surfaces with straight line path length "
                         << csf.intersection.pathLength);
      // indicate if the surface is active or not
      EX_MSG_VERBOSE(
          eCell.navigationStep,
          "surface",
          surfaceID,
          (surface->associatedDetectorElement() ? "is active" : "is passive"));
      // record the parameters as sensitive or passive depending on the surface
      ExtrapolationMode::eMode emode = surface->associatedDetectorElement()
          ? ExtrapolationMode::CollectSensitive
          : ExtrapolationMode::CollectPassive;
      // propagate to the compatible surface, return types are
      // - InProgress       : propagation to compatible surface worked
      // - Recovered        : propagation to compatible surface did not work,
      //                      leadParameters stay the same
      // - SuccessPathLimit : propagation to compatible surface
      //                      reached the path limit
      eCode = m_cfg.propagationEngine->propagate(
          eCell, *surface, pDir, {emode}, bcheck, eCell.sensitiveCurvilinear);
      CHECK_ECODE_SUCCESS_NODEST(eCell, eCode);
      // check if the propagation was successful
      if (eCode.inProgress()) {
        EX_MSG_VERBOSE(
            eCell.navigationStep,
            "layer",
            layerValue,
            "successfully hit "
                << (surface->associatedDetectorElement() ? "active" : "passive")
                << " sub structure surface.");
        // check if the surface holds material and hence needs to be processed
        if (surface->associatedMaterial()) {
          // screen output
          EX_MSG_VERBOSE(eCell.navigationStep,
                         "surface",
                         surfaceID,
                         "applying material effects.");
          // now handle the material, return codes are:
          // - SuccessMaterialLimit : material limit reached,return back
          // - InProgress           : material update done or not (depending on
          // the material description)
          eCode = m_cfg.materialEffectsEngine->handleMaterial(
              eCell, surface, pDir, fullUpdate);
          CHECK_ECODE_CONTINUE(eCell, eCode);
        }
      }
    }  // loop over test surfaces done
  }    // there are compatible surfaces

  // ----- [3] the destination situation on the layer needs to be resolved:
  // the layer is a destination layer
  // - the final propagation call is indepenent of whether sub structure was
  // resolved or not
  // - the eCell.leadParameters are at the last possible parameters
  if (sf && isFinalLayer) {
    // [B] the layer is start and destination layer but has no sub-structure
    // -> propagation to destination surface
    //  (a) the starting layer is the same layer :
    // - neither preUpdate nore postUpdate to be done, this is old-style
    // within-layer extrapolation
    // - material will be taken into account either when the layer was reached
    // from another layer
    //   or when the layer is left to another destination
    //  (b) the starting layer is not the same layer :
    // - apply the preUpdate on the parameters whein they reached the surface
    // Possible return types:
    // - SuccessDestination  : great, desintation surface hit - but post-update
    // needs to be done
    // - SuccessPathLimit    : pathlimit was reached on the way to the
    // destination surface
    eCode = m_cfg.propagationEngine->propagate(eCell,
                                               *sf,
                                               pDir,
                                               {ExtrapolationMode::Destination},
                                               false,
                                               eCell.destinationCurvilinear);
    // check for success return path limit
    CHECK_ECODE_SUCCESS_NODEST(eCell, eCode);
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "layer",
                   layerValue,
                   "attempt to hit destination surface resulted in "
                       << eCode.toString());
    // check for a potential preUpdate
    // - in case the destination surface has material and the surface was hit
    if (sf->associatedMaterial() && eCode.isSuccess()) {
      // finally do the material update, but even as this is the final call
      //  - still check for SuccessMaterialLimit
      m_cfg.materialEffectsEngine->handleMaterial(eCell, sf, pDir, preUpdate);
      // check if success was triggered through path limit reached on the way to
      // the layer
      CHECK_ECODE_SUCCESS(eCell, eCode);
    }
    // return what you have handleLayerT or extrapolateT will resolve that
    return eCode;
  }
  // return the code:
  // - if it came until here, return InProgress to not break the layer-to-layer
  return Acts::ExtrapolationCode::InProgress;
}

/// handle the failure - as configured
template <class T>
Acts::ExtrapolationCode
Acts::StaticEngine::handleReturnT(ExtrapolationCode     eCode,
                                  ExtrapolationCell<T>& eCell,
                                  const Surface*        sf,
                                  NavigationDirection   pDir,
                                  const BoundaryCheck&  bcheck) const
{
  EX_MSG_DEBUG(++eCell.navigationStep,
               "return",
               "",
               "handleReturnT with code " << eCode.toString() << " called.");
  if (eCode.isSuccessOrRecovered() || eCode.inProgress()) {
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "return",
                   "",
                   "leaving static extrapolator successfully with code "
                       << eCode.toString());
    return eCode;
  }
  EX_MSG_VERBOSE(eCell.navigationStep,
                 "return",
                 "",
                 "failure detected as "
                     << eCode.toString()
                     << " - checking fallback configuration.");
  // obviously we need a surface to exercise the fallback
  if (sf && !eCell.configurationMode(ExtrapolationMode::AvoidFallback)) {
    EX_MSG_VERBOSE(eCell.navigationStep,
                   "return",
                   "",
                   "fallback configured. Trying to "
                   "hit destination surface from last "
                   "valid parameters.");
    // check if you hit the surface, could still be stopped by PathLimit, but
    // would also count as recovered
    eCode = m_cfg.propagationEngine->propagate(eCell,
                                               *sf,
                                               pDir,
                                               {ExtrapolationMode::Propagation},
                                               bcheck,
                                               eCell.destinationCurvilinear);
  }
  // return the extrapolation code
  return eCode;
}
