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
  // @TODO check consistency
  // copy the configuration
  m_cfg = cvbConfig;
}

void
Acts::CylinderVolumeBuilder::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeBuilder::trackingVolume(
    TrackingVolumePtr  insideVolume,
    VolumeBoundsPtr    outsideBounds,
    const LayerTriple* layerTriple) const
{
  ACTS_DEBUG("Configured to build volume : " << m_cfg.volumeName);

  // the return volume
  // -----------------------------------------------------------------------------
  TrackingVolumePtr volume = nullptr;

  // now analyize the layers that are provided
  // -----------------------------------------------------
  LayerVector negativeLayers;
  LayerVector centralLayers;
  LayerVector positiveLayers;

  // - get the layers from a provided layer triple or from the layer builder
  if (layerTriple) {
    // the negative Layers
    negativeLayers = std::get<0>(*layerTriple);
    // the central Layers
    centralLayers = std::get<1>(*layerTriple);
    // the positive Layer
    positiveLayers = std::get<2>(*layerTriple);
  } else {
    // the layers are built by the layer builder
    if (m_cfg.layerBuilder) {
      // the negative Layers
      negativeLayers = m_cfg.layerBuilder->negativeLayers();
      // the central Layers
      centralLayers = m_cfg.layerBuilder->centralLayers();
      // the positive Layer
      positiveLayers = m_cfg.layerBuilder->positiveLayers();
    }
  }
  // (0) PREP WORK ------------------------------------------------
  //
  // a) inside setup
  // the volume setup for Inner
  VolumeSetup insideSetup;
  if (insideVolume) {
    // volume and inside volume
    auto ivBounds = dynamic_cast<const CylinderVolumeBounds*>(
        &insideVolume->volumeBounds());
    // set the inside values
    insideSetup.present = true;
    insideSetup.rMin    = ivBounds->innerRadius();
    insideSetup.rMax    = ivBounds->outerRadius();
    insideSetup.zMin    = insideVolume->center().z() - ivBounds->halflengthZ();
    insideSetup.zMax    = insideVolume->center().z() + ivBounds->halflengthZ();
  }
  //
  // b) outside setup
  // the volume setup for the Outside
  VolumeSetup outsideBoundSetup;
  if (outsideBounds) {
    const CylinderVolumeBounds* ocvBounds
        = dynamic_cast<const CylinderVolumeBounds*>(outsideBounds.get());
    // the cast to CylinderVolumeBounds needs to be successful
    if (ocvBounds) {
      // get values from the out bounds
      outsideBoundSetup.present = true;
      outsideBoundSetup.rMin    = ocvBounds->innerRadius();
      outsideBoundSetup.rMax    = ocvBounds->outerRadius();
      outsideBoundSetup.zMin    = -ocvBounds->halflengthZ();
      outsideBoundSetup.zMax    = ocvBounds->halflengthZ();
    }
  }

  // (A) LAYER ANALYIS ---------------------------------------------
  // analyze the layers
  VolumeSetup nVolumeSetup = analyzeLayers(negativeLayers);
  VolumeSetup cVolumeSetup = analyzeLayers(centralLayers);
  VolumeSetup pVolumeSetup = analyzeLayers(positiveLayers);
  // layer configuration
  // --------------------------------------------------------------------------
  // possbile configurations are:
  //
  // | Negagive Endcap | Barrel | Positive Endcap | -  all layers present
  //                   | Barrel |                   -  barrel present
  //                                                -  no layer present
  std::string layerConfiguration = "|";
  if (nVolumeSetup) {
    // negative layers are present
    ACTS_VERBOSE("Negative layers are present: rmin, rmax | zmin, zmax = "
                 << nVolumeSetup.toString());
    // add to the string output
    layerConfiguration += " Negative Endcap |";
  }
  if (cVolumeSetup) {
    // central layers are present
    ACTS_VERBOSE("Central layers are present:  rmin, rmax | zmin, zmax = "
                 << cVolumeSetup.toString());
    // add to the string output
    layerConfiguration += " Barrel |";
  }
  if (pVolumeSetup) {
    // positive layers are present
    ACTS_VERBOSE("Positive layers are present: rmin, rmax | zmin, zmax = "
                 << pVolumeSetup.toString());
    // add to the string output
    layerConfiguration += " Positive Endcap |";
  }
  // screen output
  ACTS_DEBUG("Layer configuration is : " << layerConfiguration);

  // (B) LAYER SETUP SYNCHRONISATION ----------------------------------
  // synchronise the layer setup
  auto wrappingSetup = synchronizeVolumeSetups(
      nVolumeSetup, cVolumeSetup, pVolumeSetup, insideSetup, outsideBoundSetup);
  if (wrappingSetup == SynchronizationError) {
    // something went wrong in the synchronisation
    ACTS_ERROR("Synchronization ERROR in layer dimensions, bailing out.")
    // return nullptr and let the upstream handle it
    return nullptr;
  }

  // (C) VOLUME CREATION ----------------------------------
  auto tvHelper = m_cfg.trackingVolumeHelper;
  // the barrel is always created
  auto barrel = tvHelper->createTrackingVolume(cVolumeSetup.layers,
                                               m_cfg.volumeMaterial,
                                               cVolumeSetup.rMin,
                                               cVolumeSetup.rMax,
                                               cVolumeSetup.zMin,
                                               cVolumeSetup.zMax,
                                               m_cfg.volumeName + "::Barrel");

  // the negative endcap is created if present
  auto nEndcap = nVolumeSetup
      ? tvHelper->createTrackingVolume(nVolumeSetup.layers,
                                       m_cfg.volumeMaterial,
                                       nVolumeSetup.rMin,
                                       nVolumeSetup.rMax,
                                       nVolumeSetup.zMin,
                                       nVolumeSetup.zMax,
                                       m_cfg.volumeName + "::NegativeEndcap")
      : nullptr;

  // the positive endcap is created
  auto pEndcap = pVolumeSetup
      ? tvHelper->createTrackingVolume(pVolumeSetup.layers,
                                       m_cfg.volumeMaterial,
                                       pVolumeSetup.rMin,
                                       pVolumeSetup.rMax,
                                       pVolumeSetup.zMin,
                                       pVolumeSetup.zMax,
                                       m_cfg.volumeName + "::PositiveEndcap")
      : nullptr;

  // (D) VOLUME WRAPPING ----------------------------------
  //
  if (wrappingSetup == NoWrapping || wrappingSetup == TripleWrapping
      || wrappingSetup == TripleWrappingGaps) {
    // we have endcap volumes
    if (nEndcap && pEndcap) {
      // a new barrel sector
      volume
          = tvHelper->createContainerTrackingVolume({nEndcap, barrel, pEndcap});
    } else  // just take the barrel as the return value
      volume = barrel;

    // the volume is now in shape for potential wrapping
    if (wrappingSetup == TripleWrappingGaps) {
      // need to create gap volume for the inside volume
      // negative gap
      auto ninGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideSetup.rMin,
          insideSetup.rMax,
          outsideBoundSetup.zMin,
          insideSetup.zMin,
          1,
          false,
          m_cfg.volumeName + "::InnerNegativeGap");
      // positive gap
      auto nipGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideSetup.rMin,
          insideSetup.rMax,
          insideSetup.zMax,
          outsideBoundSetup.zMax,
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
    if (wrappingSetup == BarrelWrappingGaps) {
      // negative gap
      auto ninGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideSetup.rMin,
          insideSetup.rMax,
          cVolumeSetup.zMin,
          insideSetup.zMin,
          1,
          false,
          m_cfg.volumeName + "::InnerNegativeGap");
      // positive gap
      auto nipGap = tvHelper->createGapTrackingVolume(
          m_cfg.volumeMaterial,
          insideSetup.rMin,
          insideSetup.rMax,
          insideSetup.zMax,
          cVolumeSetup.zMax,
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

/// synchronize the layer setups with given
Acts::CylinderVolumeBuilder::WrappingCondition
Acts::CylinderVolumeBuilder::synchronizeVolumeSetups(
    VolumeSetup&       nVolumeSetup,
    VolumeSetup&       cVolumeSetup,
    VolumeSetup&       pVolumeSetup,
    const VolumeSetup& insideSetup,
    VolumeSetup&       volumeSetup) const
{
  // for checking and estimation only
  std::vector<VolumeSetup> lsVector
      = {nVolumeSetup, cVolumeSetup, pVolumeSetup};
  // (A) estimate the inside volume if provided --------------------------
  if (insideSetup) {
    // check for potential overlap
    for (auto lSetup : lsVector)
      if (lSetup && (lSetup.overlapsInZ(insideSetup)
                     && lSetup.overlapsInR(insideSetup))) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit outside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Inside Dimensions | rmin, rmax | zmin, zmax = "
                   << insideSetup.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lSetup.toString());
        return SynchronizationError;
      }
    // everything went fine
    ACTS_VERBOSE("Inner CylinderVolumeBounds from external builder, "
                 "rmin, rmax | zmin, zmax = "
                 << insideSetup.toString());
  }

  // (B) estimate the volume dimension from outside if provided
  // -----------------
  std::string outsideMethod = "";
  if (volumeSetup) {
    for (auto lSetup : lsVector)
      if (lSetup && !volumeSetup.wraps(lSetup)) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit inside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Outside Dimensions | rmin, rmax | zmin, zmax = "
                   << volumeSetup.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lSetup.toString());
        return SynchronizationError;
      }
    // indicate the method
    outsideMethod = "provided bounds";
    // (C) outside volume provided by configuration
    // --------------------------------
  } else if (m_cfg.volumeDimension.size() > 3) {
    // cylinder volume
    volumeSetup.present = true;
    // get values from the out bounds
    volumeSetup.rMin = m_cfg.volumeDimension.at(0);
    volumeSetup.rMax = m_cfg.volumeDimension.at(1);
    volumeSetup.zMin = m_cfg.volumeDimension.at(2);
    volumeSetup.zMax = m_cfg.volumeDimension.at(3);
    // check for potential overlap
    for (auto lSetup : lsVector)
      if (lSetup && !volumeSetup.wraps(lSetup)) {
        /// give an ERROR and bail out
        ACTS_ERROR(
            "Given layer dimensions do not fit inside the provided volume "
            "bounds. Bailing out.")
        ACTS_ERROR("Outside Dimensions | rmin, rmax | zmin, zmax = "
                   << volumeSetup.toString());
        ACTS_ERROR("Layer Dimensions   | rmin, rmax | zmin, zmax = "
                   << lSetup.toString());
        return SynchronizationError;
      }
    // indicate the method
    outsideMethod = "configuration";
  } else {
    // get it from layer parsing
    for (auto lSetup : lsVector) volumeSetup.adapt(lSetup);
    // indicate the method
    outsideMethod = "layer parsing";
  }

  // screen output
  ACTS_VERBOSE("Outer bounds by " << outsideMethod
                                  << " | rmin, rmax | zmin, zmax = "
                                  << volumeSetup.toString());

  // calculate the boundary between the layers
  double zMedNeg = 0.5 * (nVolumeSetup.zMax + cVolumeSetup.zMin);
  // calculate the boundary between the layers
  double zMedPos = 0.5 * (pVolumeSetup.zMin + cVolumeSetup.zMax);

  if (nVolumeSetup && pVolumeSetup) {
    ACTS_VERBOSE("Setup | NEC || B || PEC | divisions are " << zMedNeg
                                                            << " and "
                                                            << zMedPos);
  }

  // overal dimensions are all set, now synchronize
  // always set these, there's no controversy about those
  nVolumeSetup.rMax = volumeSetup.rMax;  // n1
  pVolumeSetup.rMax = volumeSetup.rMax;  // p1
  cVolumeSetup.rMax = volumeSetup.rMax;  // c1

  nVolumeSetup.zMin = volumeSetup.zMin;  // n2
  pVolumeSetup.zMax = volumeSetup.zMax;  // p2

  auto        wCondition    = WrappingCondition::Undefined;
  std::string wConditionStr = "Undefined";

  // we need to set 4 parameters per layer setup
  // i.e. when count to 4 is done, one can bail out
  //
  // easiest case
  //
  // Case 0 - if we have no inside volume
  if (!insideSetup) {
    ACTS_VERBOSE("No inside volume provided, "
                 << " wrapping will not be neccessary.");
    // check if configured to buil to beam line
    if (m_cfg.buildToRadiusZero) volumeSetup.rMin = 0.;
    // set the minimal radius to all volume setups
    nVolumeSetup.rMin = volumeSetup.rMin;  // n3
    cVolumeSetup.rMin = volumeSetup.rMin;  // c2
    pVolumeSetup.rMin = volumeSetup.rMin;  // p3
    // 0a: no wrapping triple setup
    if (nVolumeSetup && pVolumeSetup) {
      nVolumeSetup.zMax = zMedNeg;  // n4
      cVolumeSetup.zMin = zMedNeg;  // c3
      cVolumeSetup.zMax = zMedPos;  // c4
      pVolumeSetup.zMin = zMedPos;  // p4
      // all done [ n4/c4/p4] - we can return
      wCondition = WrappingCondition::NoWrapping;
    } else {
      // 0b: only the central layer exists
      cVolumeSetup.zMin = volumeSetup.zMin;  // c3
      cVolumeSetup.zMax = volumeSetup.zMax;  // c4
      // all done [c4] - we can return
      wCondition = WrappingCondition::NoWrapping;
    }
    wConditionStr = "no Wrapping.";

  } else {
    // screen output
    ACTS_VERBOSE("Inside volume provided, "
                 << "determining wrapping condition.");
    // we have an inside volume now, one thing is for sure
    cVolumeSetup.rMin = insideSetup.rMax;  // c2

    // Case 1 - there is an inside volume
    // 1a: the inside volume is fully wrapped by the barrel
    if (cVolumeSetup.wraps(insideSetup)) {
      // the inside volume can be fully contained in the barrel
      // the outside volumes are pushed down to the inside inner radius
      nVolumeSetup.rMin = insideSetup.rMin;  // n3
      pVolumeSetup.rMin = insideSetup.rMin;  // p3
      // in that case you can set the
      nVolumeSetup.zMax = zMedNeg;  // n4
      pVolumeSetup.zMin = zMedPos;  // p4
      // set the senter volume only if pos/neg exist
      cVolumeSetup.zMin = nVolumeSetup ? zMedNeg : cVolumeSetup.zMin;  // c3
      cVolumeSetup.zMax = pVolumeSetup ? zMedPos : cVolumeSetup.zMax;  // c4
      // wrap inside the barrel,
      bool luckyBird = (cVolumeSetup.zMin == insideSetup.zMin
                        && cVolumeSetup.zMax == insideSetup.zMax);
      // but create gaps (execpt in the most lucky case)
      wCondition = luckyBird ? WrappingCondition::BarrelWrapping
                             : WrappingCondition::BarrelWrappingGaps;

      wConditionStr = luckyBird
          ? "inside barrel wrapping, barrel is already at right length."
          : "inside barrel wrapping, gaps volume to be created.";

      // 1b: the inside volume sets a comon boundary between barrel | ec
    } else if (!nVolumeSetup.overlapsInZ(insideSetup)
               && !pVolumeSetup.overlapsInZ(insideSetup)) {
      // the outside volumes are pushed down to the inside inner radius
      nVolumeSetup.rMin = insideSetup.rMin;  // n3
      pVolumeSetup.rMin = insideSetup.rMin;  // p3
      // the inside volume defines the z boundaries of the setup
      cVolumeSetup.zMin = insideSetup.zMin;  // c3
      cVolumeSetup.zMax = insideSetup.zMax;  // c4
      // and the remaining ones
      nVolumeSetup.zMax = insideSetup.zMin;  // n4
      pVolumeSetup.zMin = insideSetup.zMax;  // p4
      // wrap inside the barrel, no gaps needed, already set to right
      wCondition = WrappingCondition::BarrelWrapping;
      wConditionStr
          = "inside barrel wrapping, barrel adjusted in z to inside volume.";

      // 1c: the inside volume fits inside the triple, but needs gaps
      // 1d: the inside volume fits inside triple in R, but volume needs to be
      // extended
    } else if (nVolumeSetup && pVolumeSetup) {
      // the outside volumes are pushed down to the inside inner radius
      nVolumeSetup.rMin = insideSetup.rMin;  // n3
      pVolumeSetup.rMin = insideSetup.rMin;  // p3
      // the center extend is set
      cVolumeSetup.zMin = zMedNeg;  // c3
      cVolumeSetup.zMax = zMedPos;  // c4
      // also adjust the n/p z boundaries
      nVolumeSetup.zMax = zMedNeg;  // n4
      pVolumeSetup.zMin = zMedPos;  // p4
      // set it first to wrapping with gaps
      wCondition    = WrappingCondition::TripleWrapping;
      wConditionStr = "inside triple wrapping, gap volumes to be created.";
      // check if gaps are needed
      if (insideSetup.zMin <= nVolumeSetup.zMin
          && insideSetup.zMax >= pVolumeSetup.zMax) {
        // also need to overwrite the inner radius to close up
        nVolumeSetup.rMin = insideSetup.rMax;  // on3 - overwrites
        pVolumeSetup.rMin = insideSetup.rMax;  // op3 - overwrites
        // adjust and no wrapping needed
        nVolumeSetup.zMin = insideSetup.zMin;  // on2 - overwrites
        pVolumeSetup.zMax = insideSetup.zMax;  // op2 - overwrites
        // now overwrite it to wrappping without gaps
        wCondition = WrappingCondition::TripleWrapping;
        wConditionStr
            = "inside triple wrapping, triple adjusted in z to inside volume.";
      }
    }
  }

  // screen output after synchronization
  //
  if (insideSetup) {
    ACTS_VERBOSE("Inside volume dimensions given as rmin, rmax | zmin, zmax = "
                 << insideSetup.toString());
    ACTS_VERBOSE("Wrapping case : " << wConditionStr);
  }
  if (nVolumeSetup)
    ACTS_VERBOSE(
        "Negative volume dimensions given as rmin, rmax | zmin, zmax = "
        << nVolumeSetup.toString());
  ACTS_VERBOSE("Central volume dimensions given as  rmin, rmax | zmin, zmax = "
               << cVolumeSetup.toString());
  ACTS_VERBOSE("Positive volume dimensions given as rmin, rmax | zmin, zmax = "
               << pVolumeSetup.toString());
  ACTS_VERBOSE("Total volume dimensions given as    rmin, rmax | zmin, zmax = "
               << volumeSetup.toString());

  return wCondition;
}

// -----------------------------
Acts::VolumeSetup
Acts::CylinderVolumeBuilder::analyzeLayers(const LayerVector& lVector) const
{
  // @TODO add envelope tolerance
  //
  // return object
  VolumeSetup lSetup;
  // only if the vector is present it can actually be analyzed
  if (!lVector.empty()) {
    // we have layers
    lSetup.present = true;
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
        takeSmaller(lSetup.rMin, rMinC - m_cfg.layerEnvelopeR.first);
        takeBigger(lSetup.rMax, rMaxC + m_cfg.layerEnvelopeR.second);
        takeSmaller(lSetup.zMin, center.z() - hZ - m_cfg.layerEnvelopeZ);
        takeBigger(lSetup.zMax, center.z() + hZ + m_cfg.layerEnvelopeZ);
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
        takeSmaller(lSetup.rMin, rMinD - m_cfg.layerEnvelopeR.first);
        takeBigger(lSetup.rMax, rMaxD + m_cfg.layerEnvelopeR.second);
        takeSmaller(lSetup.zMin, zMinD - m_cfg.layerEnvelopeZ);
        takeBigger(lSetup.zMax, zMaxD + m_cfg.layerEnvelopeZ);
        //!< @TODO check for Endcap Ring setup
      }
    }
  }
  // set the layers to the layer vector
  lSetup.layers = lVector;
  // and return what you have
  return lSetup;
}
