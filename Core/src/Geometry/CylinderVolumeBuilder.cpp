// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderLayer.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/DiscLayer.hpp"
#include "Acts/Geometry/ITrackingVolumeHelper.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"

Acts::CylinderVolumeBuilder::CylinderVolumeBuilder(
    const Acts::CylinderVolumeBuilder::Config& cvbConfig,
    std::unique_ptr<const Logger> logger)
    : Acts::ITrackingVolumeBuilder(), m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(cvbConfig);
}

Acts::CylinderVolumeBuilder::~CylinderVolumeBuilder() = default;

void Acts::CylinderVolumeBuilder::setConfiguration(
    const Acts::CylinderVolumeBuilder::Config& cvbConfig) {
  // @todo check consistency
  // copy the configuration
  m_cfg = cvbConfig;
}

void Acts::CylinderVolumeBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CylinderVolumeBuilder::trackingVolume(
    const GeometryContext& gctx, TrackingVolumePtr existingVolume,
    VolumeBoundsPtr externalBounds) const {
  ACTS_DEBUG("Configured to build volume : " << m_cfg.volumeName);

  // the return volume
  // -----------------------------------------------------------------------------
  MutableTrackingVolumePtr volume = nullptr;

  // now analyize the layers that are provided
  // -----------------------------------------------------
  LayerVector negativeLayers;
  LayerVector centralLayers;
  LayerVector positiveLayers;

  // the wrapping configuration
  WrappingConfig wConfig;

  // the layers are built by the layer builder
  if (m_cfg.layerBuilder) {
    // the negative Layers
    negativeLayers = m_cfg.layerBuilder->negativeLayers(gctx);
    // the central Layers
    centralLayers = m_cfg.layerBuilder->centralLayers(gctx);
    // the positive Layer
    positiveLayers = m_cfg.layerBuilder->positiveLayers(gctx);
  }
  // (0) PREP WORK ------------------------------------------------
  //
  // a) volume config of the existing volume
  if (existingVolume) {
    // volume and existing volume
    auto existingBounds = dynamic_cast<const CylinderVolumeBounds*>(
        &existingVolume->volumeBounds());
    // set the inside values
    wConfig.existingVolumeConfig.present = true;
    wConfig.existingVolumeConfig.rMin = existingBounds->innerRadius();
    wConfig.existingVolumeConfig.rMax = existingBounds->outerRadius();
    wConfig.existingVolumeConfig.zMin =
        existingVolume->center().z() - existingBounds->halflengthZ();
    wConfig.existingVolumeConfig.zMax =
        existingVolume->center().z() + existingBounds->halflengthZ();
  }
  //
  // b) outside config
  // the volume config for the Outside
  VolumeConfig externalBoundConfig;
  if (externalBounds) {
    const CylinderVolumeBounds* ocvBounds =
        dynamic_cast<const CylinderVolumeBounds*>(externalBounds.get());
    // the cast to CylinderVolumeBounds needs to be successful
    if (ocvBounds != nullptr) {
      // get values from the out bounds
      wConfig.externalVolumeConfig.present = true;
      wConfig.externalVolumeConfig.rMin = ocvBounds->innerRadius();
      wConfig.externalVolumeConfig.rMax = ocvBounds->outerRadius();
      wConfig.externalVolumeConfig.zMin = -ocvBounds->halflengthZ();
      wConfig.externalVolumeConfig.zMax = ocvBounds->halflengthZ();
    }
  }

  // ---------------------------------------------
  // The Volume Config of the SubVolumes
  // ---------------------------------------------
  // sub volume / layer configuration (subVolumes only build of layers are
  // present)
  // --------------------------------------------------------------------------
  //
  // possbile configurations are (so far only synchronised):
  //
  // | Negative Endcap | Barrel | Positive Endcap | -  all layers present
  //                   | Barrel |                   -  barrel present
  // | Negative Endcap |        | Positive Endcap | - only endcaps present
  //                                                -  no layer present
  // Check if already given through configuration
  //
  // (A) volume configuration
  //

  // Find out with Layer analysis
  // analyze the layers
  wConfig.nVolumeConfig = analyzeLayers(gctx, negativeLayers);
  wConfig.cVolumeConfig = analyzeLayers(gctx, centralLayers);
  wConfig.pVolumeConfig = analyzeLayers(gctx, positiveLayers);

  std::string layerConfiguration = "|";
  if (wConfig.nVolumeConfig) {
    // negative layers are present
    ACTS_VERBOSE("Negative layers are present: rmin, rmax | zmin, zmax = "
                 << wConfig.nVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Negative Endcap |";
  }
  if (wConfig.cVolumeConfig) {
    // central layers are present
    ACTS_VERBOSE("Central layers are present:  rmin, rmax | zmin, zmax = "
                 << wConfig.cVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Barrel |";
  }
  if (wConfig.pVolumeConfig) {
    // positive layers are present
    ACTS_VERBOSE("Positive layers are present: rmin, rmax | zmin, zmax = "
                 << wConfig.pVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Positive Endcap |";
  }
  // screen output
  ACTS_DEBUG("Layer configuration is : " << layerConfiguration);

  // (B) LAYER Config SYNCHRONISATION ----------------------------------
  // synchronise the layer config
  ACTS_VERBOSE("Configurations after layer parsing " << '\n'
                                                     << wConfig.toString());
  // first let us arrange the new container volume
  wConfig.configureContainerVolume();
  ACTS_VERBOSE("Configuration after container synchronisation "
               << '\n'
               << wConfig.toString());
  // now let's understand the wrapping if needed
  if (wConfig.existingVolumeConfig) {
    wConfig.wrapInsertAttach();
    ACTS_VERBOSE("Configuration after wrapping, insertion, attachment "
                 << '\n'
                 << wConfig.toString());
  } else {
    // no wrapping around inner volume needed
    // however there could be central, positive & negative volume which will
    // need to be put into a container volume
    wConfig.wCondition = NoWrapping;
  }

  // (C) VOLUME CREATION ----------------------------------
  auto tvHelper = m_cfg.trackingVolumeHelper;
  // the barrel is always created
  auto barrel =
      wConfig.cVolumeConfig
          ? tvHelper->createTrackingVolume(
                gctx, wConfig.cVolumeConfig.layers, m_cfg.volumeMaterial,
                wConfig.cVolumeConfig.rMin, wConfig.cVolumeConfig.rMax,
                wConfig.cVolumeConfig.zMin, wConfig.cVolumeConfig.zMax,
                m_cfg.volumeName + "::Barrel")
          : nullptr;

  // the negative endcap is created if present
  auto nEndcap =
      wConfig.nVolumeConfig
          ? tvHelper->createTrackingVolume(
                gctx, wConfig.nVolumeConfig.layers, m_cfg.volumeMaterial,
                wConfig.nVolumeConfig.rMin, wConfig.nVolumeConfig.rMax,
                wConfig.nVolumeConfig.zMin, wConfig.nVolumeConfig.zMax,
                m_cfg.volumeName + "::NegativeEndcap")
          : nullptr;

  // the positive endcap is created
  auto pEndcap =
      wConfig.pVolumeConfig
          ? tvHelper->createTrackingVolume(
                gctx, wConfig.pVolumeConfig.layers, m_cfg.volumeMaterial,
                wConfig.pVolumeConfig.rMin, wConfig.pVolumeConfig.rMax,
                wConfig.pVolumeConfig.zMin, wConfig.pVolumeConfig.zMax,
                m_cfg.volumeName + "::PositiveEndcap")
          : nullptr;

  ACTS_DEBUG("Newly created volume(s) will be " << wConfig.wConditionScreen);
  // standalone container, full wrapping, full insertion & if no existing volume
  // is present needs a bare triple
  if (wConfig.wCondition == Wrapping || wConfig.wCondition == Inserting ||
      wConfig.wCondition == NoWrapping) {
    ACTS_VERBOSE("Combined new container is being built.");
    // stuff into the container what you have
    std::vector<std::shared_ptr<const TrackingVolume>> volumesContainer;
    if (nEndcap) {
      volumesContainer.push_back(nEndcap);
      volume = nEndcap;
    }
    if (barrel) {
      volumesContainer.push_back(barrel);
      volume = barrel;
    }
    if (pEndcap) {
      volumesContainer.push_back(pEndcap);
      volume = pEndcap;
    }
    // and low lets create the new volume
    volume =
        volumesContainer.size() > 1
            ? tvHelper->createContainerTrackingVolume(gctx, volumesContainer)
            : volume;
  } else if (wConfig.wCondition != Attaching) {
    // the new volume is the only one present
    volume = nEndcap ? nEndcap : (barrel ? barrel : pEndcap);
  }

  // prepare the gap volumes first
  TrackingVolumePtr existingVolumeCp = existingVolume;
  // check if further action is needed on existing volumes and gap volumes
  if (existingVolumeCp) {
    // check if gaps are needed
    std::vector<std::shared_ptr<const TrackingVolume>> existingContainer;
    if (wConfig.fGapVolumeConfig) {
      // create the gap volume
      auto fGap = tvHelper->createGapTrackingVolume(
          gctx, m_cfg.volumeMaterial, wConfig.fGapVolumeConfig.rMin,
          wConfig.fGapVolumeConfig.rMax, wConfig.fGapVolumeConfig.zMin,
          wConfig.fGapVolumeConfig.zMax, 1, false, m_cfg.volumeName + "::fGap");
      // push it back into the list
      existingContainer.push_back(fGap);
    }
    existingContainer.push_back(existingVolumeCp);
    if (wConfig.sGapVolumeConfig) {
      // create the gap volume
      auto sGap = tvHelper->createGapTrackingVolume(
          gctx, m_cfg.volumeMaterial, wConfig.sGapVolumeConfig.rMin,
          wConfig.sGapVolumeConfig.rMax, wConfig.sGapVolumeConfig.zMin,
          wConfig.sGapVolumeConfig.zMax, 1, false, m_cfg.volumeName + "::sGap");
      // push it back into the list
      existingContainer.push_back(sGap);
    }

    // and low lets create the new existing volume with gaps
    existingVolumeCp =
        existingContainer.size() > 1
            ? tvHelper->createContainerTrackingVolume(gctx, existingContainer)
            : existingVolumeCp;

    // for central wrapping or inserting, we need to update once more
    // clear the container
    existingContainer.clear();
    if (wConfig.wCondition == CentralWrapping) {
      existingContainer.push_back(existingVolumeCp);
      existingContainer.push_back(barrel);
    } else if (wConfig.wCondition == CentralInserting) {
      existingContainer.push_back(barrel);
      existingContainer.push_back(existingVolumeCp);
    }
    // update
    existingVolumeCp =
        !existingContainer.empty()
            ? tvHelper->createContainerTrackingVolume(gctx, existingContainer)
            : existingVolumeCp;

    std::vector<std::shared_ptr<const TrackingVolume>> totalContainer;
    // check what to do with the existing
    if (wConfig.wCondition == Attaching ||
        wConfig.wCondition == CentralWrapping ||
        wConfig.wCondition == CentralInserting) {
      if (nEndcap) {
        totalContainer.push_back(nEndcap);
      }
      totalContainer.push_back(existingVolumeCp);
      if (pEndcap) {
        totalContainer.push_back(pEndcap);
      }
    } else if (wConfig.wCondition == Inserting && volume) {
      totalContainer.push_back(volume);
      totalContainer.push_back(existingVolumeCp);
    } else if (wConfig.wCondition == Wrapping && volume) {
      totalContainer.push_back(existingVolumeCp);
      totalContainer.push_back(volume);
    } else {
      ACTS_ERROR("Misconfiguration in volume building detected.");
      return nullptr;
    }
    // now create the new container volume
    volume = tvHelper->createContainerTrackingVolume(gctx, totalContainer);
  }
  // sign the volume and return it
  volume->sign(GeometrySignature(m_cfg.volumeSignature));
  // now return what you have
  return volume;
}

// -----------------------------
Acts::VolumeConfig Acts::CylinderVolumeBuilder::analyzeLayers(
    const GeometryContext& gctx, const LayerVector& lVector) const {
  // @TODO add envelope tolerance
  //
  // return object
  VolumeConfig lConfig;
  // only if the vector is present it can actually be analyzed
  if (!lVector.empty()) {
    // we have layers
    lConfig.present = true;
    // loop over the layer
    for (auto& layer : lVector) {
      // the thickness of the layer needs to be taken into account
      double thickness = layer->thickness();
      // get the center of the layer
      const Vector3D& center = layer->surfaceRepresentation().center(gctx);
      // check if it is a cylinder layer
      const CylinderLayer* cLayer =
          dynamic_cast<const CylinderLayer*>(layer.get());
      if (cLayer != nullptr) {
        // now we have access to all the information
        double rMinC =
            cLayer->surfaceRepresentation().bounds().r() - 0.5 * thickness;
        double rMaxC =
            cLayer->surfaceRepresentation().bounds().r() + 0.5 * thickness;

        double hZ = cLayer->surfaceRepresentation().bounds().halflengthZ();
        takeSmaller(lConfig.rMin, rMinC - m_cfg.layerEnvelopeR.first);
        takeBigger(lConfig.rMax, rMaxC + m_cfg.layerEnvelopeR.second);
        takeSmaller(lConfig.zMin, center.z() - hZ - m_cfg.layerEnvelopeZ);
        takeBigger(lConfig.zMax, center.z() + hZ + m_cfg.layerEnvelopeZ);
      }
      // proceed further if it is a Disc layer
      const RadialBounds* dBounds = dynamic_cast<const RadialBounds*>(
          &(layer->surfaceRepresentation().bounds()));
      if (dBounds != nullptr) {
        // now we have access to all the information
        double rMinD = dBounds->rMin();
        double rMaxD = dBounds->rMax();
        double zMinD = center.z() - 0.5 * thickness;
        double zMaxD = center.z() + 0.5 * thickness;
        takeSmaller(lConfig.rMin, rMinD - m_cfg.layerEnvelopeR.first);
        takeBigger(lConfig.rMax, rMaxD + m_cfg.layerEnvelopeR.second);
        takeSmaller(lConfig.zMin, zMinD - m_cfg.layerEnvelopeZ);
        takeBigger(lConfig.zMax, zMaxD + m_cfg.layerEnvelopeZ);
        //!< @todo check for Endcap Ring config
      }
    }
  }
  // set the layers to the layer vector
  lConfig.layers = lVector;
  // overwrite to radius 0 if needed
  if (m_cfg.buildToRadiusZero) {
    ACTS_VERBOSE("This layer builder is configured to build to the beamline.");
    lConfig.rMin = 0.;
  }

  // and return what you have
  return lConfig;
}
