// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeBuilder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/CylinderVolumeBuilder.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Tools/ITrackingVolumeHelper.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

Acts::CylinderVolumeBuilder::CylinderVolumeBuilder(
    const Acts::CylinderVolumeBuilder::Config& cvbConfig,
    std::unique_ptr<Logger>                    logger)
  : Acts::ITrackingVolumeBuilder(), m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(cvbConfig);
}

Acts::CylinderVolumeBuilder::~CylinderVolumeBuilder()
{
}

void
Acts::CylinderVolumeBuilder::setConfiguration(
    const Acts::CylinderVolumeBuilder::Config& cvbConfig)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = cvbConfig;
}

void
Acts::CylinderVolumeBuilder::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

std::shared_ptr<Acts::TrackingVolume>
Acts::CylinderVolumeBuilder::trackingVolume(TrackingVolumePtr insideVolume,
                                            VolumeBoundsPtr outsideBounds) const
{
  ACTS_DEBUG("Configured to build volume : " << m_cfg.volumeName);

  // the return volume
  // -----------------------------------------------------------------------------
  MutableTrackingVolumePtr volume = nullptr;

  // now analyize the layers that are provided
  // -----------------------------------------------------
  LayerVector negativeLayers;
  LayerVector centralLayers;
  LayerVector positiveLayers;

  // the layers are built by the layer builder
  if (m_cfg.layerBuilder) {
    // the negative Layers
    negativeLayers = m_cfg.layerBuilder->negativeLayers();
    // the central Layers
    centralLayers = m_cfg.layerBuilder->centralLayers();
    // the positive Layer
    positiveLayers = m_cfg.layerBuilder->positiveLayers();
  }
  // (0) PREP WORK ------------------------------------------------
  //
  // a) inside config
  // the volume config for Inner
  VolumeConfig insideConfig;
  if (insideVolume) {
    // volume and inside volume
    auto ivBounds = dynamic_cast<const CylinderVolumeBounds*>(
        &insideVolume->volumeBounds());
    // set the inside values
    insideConfig.present = true;
    insideConfig.rMin    = ivBounds->innerRadius();
    insideConfig.rMax    = ivBounds->outerRadius();
    insideConfig.zMin    = insideVolume->center().z() - ivBounds->halflengthZ();
    insideConfig.zMax    = insideVolume->center().z() + ivBounds->halflengthZ();
  }
  //
  // b) outside config
  // the volume config for the Outside
  VolumeConfig outsideBoundConfig;
  if (outsideBounds) {
    const CylinderVolumeBounds* ocvBounds
        = dynamic_cast<const CylinderVolumeBounds*>(outsideBounds.get());
    // the cast to CylinderVolumeBounds needs to be successful
    if (ocvBounds) {
      // get values from the out bounds
      outsideBoundConfig.present = true;
      outsideBoundConfig.rMin    = ocvBounds->innerRadius();
      outsideBoundConfig.rMax    = ocvBounds->outerRadius();
      outsideBoundConfig.zMin    = -ocvBounds->halflengthZ();
      outsideBoundConfig.zMax    = ocvBounds->halflengthZ();
    }
  }
  // ---------------------------------------------
  // The Volume Config of the SubVolumes
  // ---------------------------------------------
  // sub volume / layer configuration (subVolumes only build of layers are
  // present)
  // --------------------------------------------------------------------------
  // possbile configurations are:
  //
  // | Negative Endcap | Barrel | Positive Endcap | -  all layers present
  //                   | Barrel |                   -  barrel present
  //                                                -  no layer present
  VolumeConfig nVolumeConfig;
  VolumeConfig cVolumeConfig;
  VolumeConfig pVolumeConfig;
  // Check if already given
  if (m_cfg.subVolumeConfig) {
    if (!negativeLayers.empty()) {
      if (m_cfg.subVolumeConfig.zBoundaries.size() < 4) {
        ACTS_ERROR("Only " << m_cfg.subVolumeConfig.zBoundaries.size()
                           << " zBoundaries are given for the subVolumeConfig, "
                              "but negative Layers are present and 4 "
                              "zBoundaries need to be given in this case - "
                              "Please check the configuration!"
                              "The negative Volume will not be built now.");
        nVolumeConfig.present = false;
      } else {
        // we have layers
        nVolumeConfig.present = true;
        nVolumeConfig.rMin    = m_cfg.subVolumeConfig.rMin;
        nVolumeConfig.rMax    = m_cfg.subVolumeConfig.rMax;
        nVolumeConfig.zMin    = m_cfg.subVolumeConfig.zBoundaries.at(0);
        nVolumeConfig.zMax    = m_cfg.subVolumeConfig.zBoundaries.at(1);
        nVolumeConfig.layers  = negativeLayers;
      }
    }
    if (!centralLayers.empty()) {
      // we have layers
      cVolumeConfig.present = true;
      cVolumeConfig.rMin    = m_cfg.subVolumeConfig.rMin;
      cVolumeConfig.rMax    = m_cfg.subVolumeConfig.rMax;
      cVolumeConfig.zMin    = (m_cfg.subVolumeConfig.zBoundaries.size() < 4)
          ? m_cfg.subVolumeConfig.zBoundaries.at(0)
          : m_cfg.subVolumeConfig.zBoundaries.at(1);
      cVolumeConfig.zMax = (m_cfg.subVolumeConfig.zBoundaries.size() < 4)
          ? m_cfg.subVolumeConfig.zBoundaries.at(1)
          : m_cfg.subVolumeConfig.zBoundaries.at(2);
      cVolumeConfig.layers = centralLayers;
      cVolumeConfig.layers = centralLayers;
    }
    if (!positiveLayers.empty()) {
      if (m_cfg.subVolumeConfig.zBoundaries.size() < 4) {
        ACTS_ERROR("Only " << m_cfg.subVolumeConfig.zBoundaries.size()
                           << " zBoundaries are given for the subVolumeConfig, "
                              "but positive Layers are present and 4 "
                              "zBoundaries need to be given in this case - "
                              "Please check the configuration!"
                              "The positive Volume will not be built now.");
        pVolumeConfig.present = false;
      } else {
        // we have layers
        pVolumeConfig.present = true;
        pVolumeConfig.rMin    = m_cfg.subVolumeConfig.rMin;
        pVolumeConfig.rMax    = m_cfg.subVolumeConfig.rMax;
        pVolumeConfig.zMin    = m_cfg.subVolumeConfig.zBoundaries.at(2);
        pVolumeConfig.zMax    = m_cfg.subVolumeConfig.zBoundaries.at(3);
        pVolumeConfig.layers  = positiveLayers;
      }
    }
  } else {
    // Find out with Layer analysis
    // (A) LAYER ANALYIS ---------------------------------------------
    // analyze the layers
    nVolumeConfig = analyzeLayers(negativeLayers);
    cVolumeConfig = analyzeLayers(centralLayers);
    pVolumeConfig = analyzeLayers(positiveLayers);
  }
  std::string layerConfiguration = "|";
  if (nVolumeConfig) {
    // negative layers are present
    ACTS_VERBOSE("Negative layers are present: rmin, rmax | zmin, zmax = "
                 << nVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Negative Endcap |";
  }
  if (cVolumeConfig) {
    // central layers are present
    ACTS_VERBOSE("Central layers are present:  rmin, rmax | zmin, zmax = "
                 << cVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Barrel |";
  }
  if (pVolumeConfig) {
    // positive layers are present
    ACTS_VERBOSE("Positive layers are present: rmin, rmax | zmin, zmax = "
                 << pVolumeConfig.toString());
    // add to the string output
    layerConfiguration += " Positive Endcap |";
  }
  // screen output
  ACTS_DEBUG("Layer configuration is : " << layerConfiguration);

  // (B) LAYER Config SYNCHRONISATION ----------------------------------
  // synchronise the layer config
  auto wrappingConfig = synchronizeVolumeConfigs(nVolumeConfig,
                                                 cVolumeConfig,
                                                 pVolumeConfig,
                                                 insideConfig,
                                                 outsideBoundConfig);
  if (wrappingConfig == SynchronizationError) {
    // something went wrong in the synchronisation
    ACTS_ERROR("Synchronization ERROR in layer dimensions, bailing out.")
    // return nullptr and let the upstream handle it
    return nullptr;
  }

  // (C) VOLUME CREATION ----------------------------------
  auto tvHelper = m_cfg.trackingVolumeHelper;
  // the barrel is always created
  auto barrel = tvHelper->createTrackingVolume(cVolumeConfig.layers,
                                               m_cfg.volumeMaterial,
                                               cVolumeConfig.rMin,
                                               cVolumeConfig.rMax,
                                               cVolumeConfig.zMin,
                                               cVolumeConfig.zMax,
                                               m_cfg.volumeName + "::Barrel");

  // the negative endcap is created if present
  auto nEndcap = nVolumeConfig
      ? tvHelper->createTrackingVolume(nVolumeConfig.layers,
                                       m_cfg.volumeMaterial,
                                       nVolumeConfig.rMin,
                                       nVolumeConfig.rMax,
                                       nVolumeConfig.zMin,
                                       nVolumeConfig.zMax,
                                       m_cfg.volumeName + "::NegativeEndcap")
      : nullptr;

  // the positive endcap is created
  auto pEndcap = pVolumeConfig
      ? tvHelper->createTrackingVolume(pVolumeConfig.layers,
                                       m_cfg.volumeMaterial,
                                       pVolumeConfig.rMin,
                                       pVolumeConfig.rMax,
                                       pVolumeConfig.zMin,
                                       pVolumeConfig.zMax,
                                       m_cfg.volumeName + "::PositiveEndcap")
      : nullptr;

  // (D) VOLUME WRAPPING ----------------------------------
  //
  if (wrappingConfig == NoWrapping || wrappingConfig == TripleWrapping
      || wrappingConfig == TripleWrappingGaps) {
    // we have endcap volumes
    if (nEndcap && pEndcap) {
      // a new barrel sector
      volume
          = tvHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
    } else  // just take the barrel as the return value
      volume = barrel;

    // the volume is now in shape for potential wrapping
    if (wrappingConfig == TripleWrappingGaps) {
      // need to create gap volume for the inside volume
      // negative gap
      auto ninGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideConfig.rMin,
          insideConfig.rMax,
          outsideBoundConfig.zMin,
          insideConfig.zMin,
          1,
          false,
          m_cfg.volumeName + "::InnerNegativeGap");
      // positive gap
      auto nipGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideConfig.rMin,
          insideConfig.rMax,
          insideConfig.zMax,
          outsideBoundConfig.zMax,
          1,
          false,
          m_cfg.volumeName + "::InnerPositiveGap");
      // update the inside volume to the right dimensions
      insideVolume = tvHelper->createContainerTrackingVolume(
          {ninGap, insideVolume, nipGap});
    }
    // if there is an inside volume, then wrap - otherwise : done
    volume = insideVolume
        ? tvHelper->createContainerTrackingVolume({insideVolume, volume})
        : volume;
    // done
  } else {
    // the inside volume may have to be adapted to the barrel length
    if (wrappingConfig == BarrelWrappingGaps) {
      // negative gap
      auto ninGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideConfig.rMin,
          insideConfig.rMax,
          cVolumeConfig.zMin,
          insideConfig.zMin,
          1,
          false,
          m_cfg.volumeName + "::InnerNegativeGap");
      // positive gap
      auto nipGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideConfig.rMin,
          insideConfig.rMax,
          insideConfig.zMax,
          cVolumeConfig.zMax,
          1,
          false,
          m_cfg.volumeName + "::InnerPositiveGap");
      // update the inside volume to the right dimensions
      insideVolume = tvHelper->createContainerTrackingVolume(
          {ninGap, insideVolume, nipGap});
    }
    // now wrap the insideVolume into the Barrel
    barrel = tvHelper->createContainerTrackingVolume({insideVolume, barrel});

    // and now update the volume with the endcaps (if exist)
    if (nEndcap && pEndcap) {
      // a new barrel sector
      volume
          = tvHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
    } else  // just take the barrel as the return value
      volume = barrel;
  }

  // sign the volume and return it
  volume->sign(GeometrySignature(m_cfg.volumeSignature));
  // now return what you have
  return volume;
}

/// synchronize the layer configs with given
Acts::CylinderVolumeBuilder::WrappingCondition
Acts::CylinderVolumeBuilder::synchronizeVolumeConfigs(
    VolumeConfig&       nVolumeConfig,
    VolumeConfig&       cVolumeConfig,
    VolumeConfig&       pVolumeConfig,
    const VolumeConfig& insideConfig,
    VolumeConfig&       volumeConfig) const
{
  // for checking and estimation only
  std::vector<VolumeConfig> lsVector
      = {nVolumeConfig, cVolumeConfig, pVolumeConfig};
  // (A) estimate the inside volume if provided --------------------------
  if (insideConfig) {
    // check for potential overlap
    for (auto lConfig : lsVector)
      if (lConfig && !lConfig.wraps(insideConfig)) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit outside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Inside Dimensions | rmin, rmax | zmin, zmax = "
                   << insideConfig.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lConfig.toString());
        return SynchronizationError;
      }
    // everything went fine
    ACTS_VERBOSE("Inner CylinderVolumeBounds from external builder, "
                 "rmin, rmax | zmin, zmax = "
                 << insideConfig.toString());
  }

  // (B) estimate the volume dimension from outside if provided
  // -----------------
  std::string outsideMethod = "";
  if (volumeConfig) {
    for (auto lConfig : lsVector)
      if (lConfig && !lConfig.containes(volumeConfig)) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit inside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Outside Dimensions | rmin, rmax | zmin, zmax = "
                   << volumeConfig.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lConfig.toString());
        return SynchronizationError;
      }
    // indicate the method
    outsideMethod = "provided bounds";
    // (C) outside volume provided by configuration
    // --------------------------------
  } else if (m_cfg.volumeDimension.size() > 3) {
    // cylinder volume
    volumeConfig.present = true;
    // get values from the out bounds
    volumeConfig.rMin = m_cfg.volumeDimension.at(0);
    volumeConfig.rMax = m_cfg.volumeDimension.at(1);
    volumeConfig.zMin = m_cfg.volumeDimension.at(2);
    volumeConfig.zMax = m_cfg.volumeDimension.at(3);
    // check for potential overlap
    for (auto lConfig : lsVector)
      if (lConfig && !lConfig.containes(volumeConfig)) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit inside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Outside Dimensions | rmin, rmax | zmin, zmax = "
                   << volumeConfig.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lConfig.toString());
        return SynchronizationError;
      }
    // indicate the method
    outsideMethod = "configuration";
  } else if (m_cfg.subVolumeConfig) {
    volumeConfig.present = true;
    // get the bounds from teh sub volumes
    volumeConfig.rMin = m_cfg.subVolumeConfig.rMin;
    volumeConfig.rMax = m_cfg.subVolumeConfig.rMax;
    volumeConfig.zMin = m_cfg.subVolumeConfig.zBoundaries.front();
    volumeConfig.zMax = m_cfg.subVolumeConfig.zBoundaries.back();
  } else {
    // get it from layer parsing
    for (auto lConfig : lsVector) volumeConfig.adapt(lConfig);
    // indicate the method
    outsideMethod = "layer parsing";
  }

  // screen output
  ACTS_VERBOSE("Outer bounds by " << outsideMethod
                                  << " | rmin, rmax | zmin, zmax = "
                                  << volumeConfig.toString());

  // calculate the boundary between the layers
  double zMedNeg = 0.5 * (nVolumeConfig.zMax + cVolumeConfig.zMin);
  // calculate the boundary between the layers
  double zMedPos = 0.5 * (pVolumeConfig.zMin + cVolumeConfig.zMax);

  if (nVolumeConfig && pVolumeConfig) {
    ACTS_VERBOSE("Config | NEC || B || PEC | divisions are " << zMedNeg
                                                             << " and "
                                                             << zMedPos);
  }

  // overal dimensions are all set, now synchronize
  // always set these, there's no controversy about those
  nVolumeConfig.rMax = volumeConfig.rMax;  // n1
  pVolumeConfig.rMax = volumeConfig.rMax;  // p1
  cVolumeConfig.rMax = volumeConfig.rMax;  // c1

  nVolumeConfig.zMin = volumeConfig.zMin;  // n2
  pVolumeConfig.zMax = volumeConfig.zMax;  // p2

  auto        wCondition    = WrappingCondition::Undefined;
  std::string wConditionStr = "Undefined";

  // we need to set 4 parameters per layer Config
  // i.e. when count to 4 is done, one can bail out
  //
  // easiest case
  //
  // Case 0 - if we have no inside volume
  if (!insideConfig) {
    ACTS_VERBOSE("No inside volume provided, "
                 << " wrapping will not be neccessary.");
    // check if configured to buil to beam line
    if (m_cfg.buildToRadiusZero) volumeConfig.rMin = 0.;
    // set the minimal radius to all volume configs
    nVolumeConfig.rMin = volumeConfig.rMin;  // n3
    cVolumeConfig.rMin = volumeConfig.rMin;  // c2
    pVolumeConfig.rMin = volumeConfig.rMin;  // p3
    // 0a: no wrapping triple Config
    if (nVolumeConfig && pVolumeConfig) {
      nVolumeConfig.zMax = zMedNeg;  // n4
      cVolumeConfig.zMin = zMedNeg;  // c3
      cVolumeConfig.zMax = zMedPos;  // c4
      pVolumeConfig.zMin = zMedPos;  // p4
      // all done [ n4/c4/p4] - we can return
      wCondition = WrappingCondition::NoWrapping;
    } else {
      // 0b: only the central layer exists
      cVolumeConfig.zMin = volumeConfig.zMin;  // c3
      cVolumeConfig.zMax = volumeConfig.zMax;  // c4
      // all done [c4] - we can return
      wCondition = WrappingCondition::NoWrapping;
    }
    wConditionStr = "no Wrapping.";

  } else {
    // screen output
    ACTS_VERBOSE("Inside volume provided, "
                 << "determining wrapping condition.");
    // we have an inside volume now, one thing is for sure
    cVolumeConfig.rMin = insideConfig.rMax;  // c2

    // Case 1 - there is an inside volume
    // 1a: the inside volume is fully wrapped by the barrel
    if (cVolumeConfig.containes(insideConfig)) {
      // the inside volume can be fully contained in the barrel
      // the outside volumes are pushed down to the inside inner radius
      nVolumeConfig.rMin = insideConfig.rMin;  // n3
      pVolumeConfig.rMin = insideConfig.rMin;  // p3
      // in that case you can set the
      nVolumeConfig.zMax = zMedNeg;  // n4
      pVolumeConfig.zMin = zMedPos;  // p4
      // set the senter volume only if pos/neg exist
      cVolumeConfig.zMin = nVolumeConfig ? zMedNeg : cVolumeConfig.zMin;  // c3
      cVolumeConfig.zMax = pVolumeConfig ? zMedPos : cVolumeConfig.zMax;  // c4
      // wrap inside the barrel,
      bool luckyBird = (cVolumeConfig.zMin == insideConfig.zMin
                        && cVolumeConfig.zMax == insideConfig.zMax);
      // but create gaps (execpt in the most lucky case)
      wCondition = luckyBird ? WrappingCondition::BarrelWrapping
                             : WrappingCondition::BarrelWrappingGaps;

      wConditionStr = luckyBird
          ? "inside barrel wrapping, barrel is already at right length."
          : "inside barrel wrapping, gaps volume to be created.";

      // 1b: the inside volume sets a comon boundary between barrel | ec
    } else if (!nVolumeConfig.overlapsInZ(insideConfig)
               && !pVolumeConfig.overlapsInZ(insideConfig)) {
      // the outside volumes are pushed down to the inside inner radius
      nVolumeConfig.rMin = insideConfig.rMin;  // n3
      pVolumeConfig.rMin = insideConfig.rMin;  // p3
      // the inside volume defines the z boundaries of the config
      cVolumeConfig.zMin = insideConfig.zMin;  // c3
      cVolumeConfig.zMax = insideConfig.zMax;  // c4
      // and the remaining ones
      nVolumeConfig.zMax = insideConfig.zMin;  // n4
      pVolumeConfig.zMin = insideConfig.zMax;  // p4
      // wrap inside the barrel, no gaps needed, already set to right
      wCondition = WrappingCondition::BarrelWrapping;
      wConditionStr
          = "inside barrel wrapping, barrel adjusted in z to inside volume.";

      // 1c: the inside volume fits inside the triple, but needs gaps
      // 1d: the inside volume fits inside triple in R, but volume needs to be
      // extended
    } else if (nVolumeConfig && pVolumeConfig) {
      // the outside volumes are pushed down to the inside inner radius
      nVolumeConfig.rMin = insideConfig.rMin;  // n3
      pVolumeConfig.rMin = insideConfig.rMin;  // p3
      // the center extend is set
      cVolumeConfig.zMin = zMedNeg;  // c3
      cVolumeConfig.zMax = zMedPos;  // c4
      // also adjust the n/p z boundaries
      nVolumeConfig.zMax = zMedNeg;  // n4
      pVolumeConfig.zMin = zMedPos;  // p4
      // set it first to wrapping with gaps
      wCondition    = WrappingCondition::TripleWrapping;
      wConditionStr = "inside triple wrapping, gap volumes to be created.";
      // check if gaps are needed
      if (insideConfig.zMin <= nVolumeConfig.zMin
          && insideConfig.zMax >= pVolumeConfig.zMax) {
        // also need to overwrite the inner radius to close up
        nVolumeConfig.rMin = insideConfig.rMax;  // on3 - overwrites
        pVolumeConfig.rMin = insideConfig.rMax;  // op3 - overwrites
        // adjust and no wrapping needed
        nVolumeConfig.zMin = insideConfig.zMin;  // on2 - overwrites
        pVolumeConfig.zMax = insideConfig.zMax;  // op2 - overwrites
        // now overwrite it to wrappping without gaps
        wCondition = WrappingCondition::TripleWrapping;
        wConditionStr
            = "inside triple wrapping, triple adjusted in z to inside volume.";
      }
    }
  }

  // screen output after synchronization
  //
  if (insideConfig) {
    ACTS_VERBOSE("Inside volume dimensions given as rmin, rmax | zmin, zmax = "
                 << insideConfig.toString());
    ACTS_VERBOSE("Wrapping case : " << wConditionStr);
  }
  if (nVolumeConfig)
    ACTS_VERBOSE(
        "Negative volume dimensions given as rmin, rmax | zmin, zmax = "
        << nVolumeConfig.toString());
  ACTS_VERBOSE("Central volume dimensions given as  rmin, rmax | zmin, zmax = "
               << cVolumeConfig.toString());
  ACTS_VERBOSE("Positive volume dimensions given as rmin, rmax | zmin, zmax = "
               << pVolumeConfig.toString());
  ACTS_VERBOSE("Total volume dimensions given as    rmin, rmax | zmin, zmax = "
               << volumeConfig.toString());

  return wCondition;
}

// -----------------------------
Acts::VolumeConfig
Acts::CylinderVolumeBuilder::analyzeLayers(const LayerVector& lVector) const
{
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
      const Vector3D& center = layer->surfaceRepresentation().center();
      // check if it is a cylinder layer
      const CylinderLayer* cLayer
          = dynamic_cast<const CylinderLayer*>(layer.get());
      if (cLayer) {
        // now we have access to all the information
        double rMinC
            = cLayer->surfaceRepresentation().bounds().r() - 0.5 * thickness;
        double rMaxC
            = cLayer->surfaceRepresentation().bounds().r() + 0.5 * thickness;

        double hZ = cLayer->surfaceRepresentation().bounds().halflengthZ();
        takeSmaller(lConfig.rMin, rMinC - m_cfg.layerEnvelopeR.first);
        takeBigger(lConfig.rMax, rMaxC + m_cfg.layerEnvelopeR.second);
        takeSmaller(lConfig.zMin, center.z() - hZ - m_cfg.layerEnvelopeZ);
        takeBigger(lConfig.zMax, center.z() + hZ + m_cfg.layerEnvelopeZ);
      }
      // proceed further if it is a Disc layer
      const RadialBounds* dBounds = dynamic_cast<const RadialBounds*>(
          &(layer->surfaceRepresentation().bounds()));
      if (dBounds) {
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
  // and return what you have
  return lConfig;
}
