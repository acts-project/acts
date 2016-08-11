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
#include "ACTS/Utilities/MsgMacros.hpp"
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
    TrackingVolumePtr        insideVolume,
    VolumeBoundsPtr          outsideBounds,
    const Acts::LayerTriple* layerTriple,
    const VolumeTriple*      volumeTriple) const
{
  MSG_DEBUG("Configured to build volume : " << m_cfg.volumeName);

  // the return volume
  // -----------------------------------------------------------------------------
  std::shared_ptr<const TrackingVolume> volume = nullptr;
  // used throughout
  TrackingVolumePtr nEndcap = nullptr;
  TrackingVolumePtr barrel  = nullptr;
  TrackingVolumePtr pEndcap = nullptr;

  // get the full extend
  double volumeRmin = 10e10;
  double volumeRmax = -10e10;
  double volumeZmax = 0.;
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
    if (m_cfg.layerBuilder) {
      // the negative Layers
      negativeLayers = m_cfg.layerBuilder->negativeLayers();
      // the central Layers
      centralLayers = m_cfg.layerBuilder->centralLayers();
      // the positive Layer
      positiveLayers = m_cfg.layerBuilder->positiveLayers();
    }
  }
  // analyze the layers
  LayerSetup nLayerSetup = analyzeLayerSetup(negativeLayers);
  LayerSetup cLayerSetup = analyzeLayerSetup(centralLayers);
  LayerSetup pLayerSetup = analyzeLayerSetup(positiveLayers);
  // layer configuration
  // --------------------------------------------------------------------------
  // dimensions
  double layerRmin = 10e10;
  double layerRmax = 0.;
  double layerZmax = 0.;
  // possbile configurations are:
  //
  // 111 - all     layers present
  // 010 - central layers present
  //  0  - no layers present
  int layerConfiguration = 0;
  if (nLayerSetup) {
    // negative layers are present
    ACTS_DEBUG("Negative layers are present with r(min,max) = ( "
               << nLayerSetup.rBoundaries.first
               << ", "
               << nLayerSetup.rBoundaries.second
               << " )");
    ACTS_DEBUG("                                 z(min,max) = ( "
               << nLayerSetup.zBoundaries.first
               << ", "
               << nLayerSetup.zBoundaries.second
               << " )");

    takeSmaller(layerRmin, nLayerSetup.rBoundaries.first);
    takeBigger(layerRmax, nLayerSetup.rBoundaries.second);
    takeBigger(layerZmax, fabs(nLayerSetup.zBoundaries.first));
    // set the 100-digit for n present
    layerConfiguration += 100;
  }
  if (cLayerSetup) {
    // central layers are present
    ACTS_DEBUG("Central  layers are present with r(min,max) = ( "
               << cLayerSetup.rBoundaries.first
               << ", "
               << cLayerSetup.rBoundaries.second
               << " )");
    ACTS_DEBUG("                                 z(min,max) = ( "
               << cLayerSetup.zBoundaries.first
               << ", "
               << cLayerSetup.zBoundaries.second
               << " )");

    takeSmaller(layerRmin, cLayerSetup.rBoundaries.first);
    takeBigger(layerRmax, cLayerSetup.rBoundaries.second);
    takeBigger(layerZmax, fabs(cLayerSetup.zBoundaries.first));
    takeBigger(layerZmax, cLayerSetup.zBoundaries.second);
    // set the 10-digit for c present
    layerConfiguration += 10;
  }
  if (pLayerSetup) {
    // positive layers are present
    ACTS_DEBUG("Positive layers are present with r(min,max) = ( "
               << pLayerSetup.rBoundaries.first
               << ", "
               << pLayerSetup.rBoundaries.second
               << " )");
    ACTS_DEBUG("                                 z(min,max) = ( "
               << pLayerSetup.zBoundaries.first
               << ", "
               << pLayerSetup.zBoundaries.second
               << " )");

    takeSmaller(layerRmin, pLayerSetup.rBoundaries.first);
    takeBigger(layerRmax, pLayerSetup.rBoundaries.second);
    takeBigger(layerZmax, pLayerSetup.zBoundaries.second);
    // set the 1-digit for p present
    layerConfiguration += 1;
  }

  ACTS_DEBUG("Layer configuration estimated with " << layerConfiguration);

  // the inside volume dimensions
  // ------------------------------------------------------------------
  double insideVolumeRmin = 0.;
  double insideVolumeRmax = 0.;
  double insideVolumeZmax = 0.;
  if (insideVolume) {
    // cast to cylinder volume
    const CylinderVolumeBounds* icvBounds
        = dynamic_cast<const CylinderVolumeBounds*>(
            &(insideVolume->volumeBounds()));
    // cylindrical volume bounds are there
    if (icvBounds) {
      // the outer radius of the inner volume
      insideVolumeRmin = icvBounds->innerRadius();
      insideVolumeRmax = icvBounds->outerRadius();
      insideVolumeZmax = insideVolume->center().z() + icvBounds->halflengthZ();
      MSG_VERBOSE("Inner CylinderVolumeBounds provided from external builder, "
                  "rMin/rMax/zMax = "
                  << insideVolumeRmin
                  << ", "
                  << insideVolumeRmax
                  << ", "
                  << insideVolumeZmax);
    } else {
      // we need to bail out, the given volume is not cylindrical
      ACTS_ERROR("Given volume to wrap was not cylindrical. Bailing out.");
      // cleanup teh memory
      negativeLayers.clear();
      centralLayers.clear();
      positiveLayers.clear();
      // return a null pointer, upstream builder will have to understand this
      return nullptr;
    }
  }
  // -------------------- outside boundary conditions
  // --------------------------------------------------
  // check if we have outsideBounds
  if (outsideBounds) {
    const CylinderVolumeBounds* ocvBounds
        = dynamic_cast<const CylinderVolumeBounds*>(outsideBounds.get());
    // the cast to CylinderVolumeBounds needs to be successful
    if (ocvBounds) {
      // get values from the out bounds
      volumeRmin = ocvBounds->innerRadius();
      volumeRmax = ocvBounds->outerRadius();
      volumeZmax = ocvBounds->halflengthZ();
      ACTS_VERBOSE("Outer CylinderVolumeBounds provided from external builder, "
                   "rMin/rMax/zMax = "
                   << volumeRmin
                   << ", "
                   << volumeRmax
                   << ", "
                   << volumeZmax);
    } else {
      ACTS_ERROR("Non-cylindrical bounds given to the CylinderVolumeBuilder. "
                 "Bailing out.");
      // cleanup teh memory
      negativeLayers.clear();
      centralLayers.clear();
      positiveLayers.clear();
      // return a null pointer, upstream builder will have to understand this
      return nullptr;
    }
    // check if the outside bounds cover all the layers
    if (layerConfiguration && (volumeRmin > layerRmin || volumeRmax < layerRmax
                               || volumeZmax < layerZmax)) {
      ACTS_ERROR("Given layer dimensions do not fit inside the provided volume "
                 "bounds. Bailing out."
                 << " volumeRmin: "
                 << volumeRmin
                 << " volumeRmax: "
                 << volumeRmax
                 << " layerRmin: "
                 << layerRmin
                 << " layerRmax: "
                 << layerRmax
                 << " volumeZmax: "
                 << volumeZmax
                 << " layerZmax: "
                 << layerZmax);
      // cleanup teh memory
      negativeLayers.clear();
      centralLayers.clear();
      positiveLayers.clear();
      // return a null pointer, upstream builder will have to understand this
      return nullptr;
    }
  } else if (m_cfg.volumeDimension.size() > 2) {
    // cylinder volume
    // get values from the out bounds
    volumeRmin = m_cfg.volumeDimension.at(0);
    volumeRmax = m_cfg.volumeDimension.at(1);
    volumeZmax = m_cfg.volumeDimension.at(2);
    ACTS_VERBOSE("Outer CylinderVolumeBounds provided by configuration, "
                 "rMin/rMax/zMax = "
                 << volumeRmin
                 << ", "
                 << volumeRmax
                 << ", "
                 << volumeZmax);
  } else {
    // outside dimensions will have to be determined by the layer dimensions
    volumeRmin = m_cfg.volumeToBeamPipe ? 0. : layerRmin - m_cfg.layerEnvelopeR;
    volumeRmax = layerRmax + m_cfg.layerEnvelopeR;
    volumeZmax = layerZmax + m_cfg.layerEnvelopeZ;
    // from setup
    ACTS_VERBOSE("Outer CylinderVolumeBounds estimated from layer setup, "
                 "rMin/rMax/zMax = "
                 << volumeRmin
                 << ", "
                 << volumeRmax
                 << ", "
                 << volumeZmax);
  }
  // -------------------- analyse the layer setups
  // --------------------------------------------------
  TrackingVolumePtr negativeSector = nullptr;
  TrackingVolumePtr centralSector  = nullptr;
  TrackingVolumePtr positiveSector = nullptr;

  // wrapping condition
  // 0 - no wrapping
  // 1 - wrap central barrel plus endcap volumes around inside volume
  //   - (a) gap volumes will have to be created to extend to potential z extend
  //   (if required)
  // 2 - wrap full setup around inside volume (fitting)
  //   - (a) barrel without endcap volumes
  //   - (b) endcaps are present and their position in z is around the inside
  //   volume
  //   - (c) gap volumes will have to be created to extend to potential z extend
  //   (if required)
  int wrappingCondition = 0;
  // check if layers are present
  if (layerConfiguration) {
    // screen output
    ACTS_DEBUG("Building Volume from layer configuration.");
    // barrel configuration
    double barrelRmin = 0.;
    double barrelRmax = 0.;
    double barrelZmax = 0.;
    // endcap configuration
    double endcapRmin = 0.;
    double endcapRmax = 0.;
    double endcapZmin = 0.;
    double endcapZmax = 0.;
    // if the containing volumes are given, get the boundaries of them
    if (volumeTriple) {
      VolumePtr nEndcapVolume = std::get<0>(*volumeTriple);
      VolumePtr barrelVolume  = std::get<1>(*volumeTriple);
      VolumePtr endcapVolume  = std::get<2>(*volumeTriple);
      if (barrelVolume) {
        const CylinderVolumeBounds* barrelBounds
            = dynamic_cast<const CylinderVolumeBounds*>(
                &(barrelVolume->volumeBounds()));
        barrelRmin = barrelBounds->innerRadius();
        barrelRmax = barrelBounds->outerRadius();
        barrelZmax = barrelVolume->center().z() + barrelBounds->halflengthZ();

        ACTS_VERBOSE(
            "Outer Barrel bounds provided by configuration, rMin/rMax/zMax = "
            << barrelRmin
            << ", "
            << barrelRmax
            << ", "
            << barrelZmax);
      } else
        ACTS_ERROR("No Barrel volume given for current hierarchy!");

      // check if end cap volumes are provided
      if (endcapVolume) {
        const CylinderVolumeBounds* endcapBounds
            = dynamic_cast<const CylinderVolumeBounds*>(
                &(endcapVolume->volumeBounds()));
        endcapRmin = endcapBounds->innerRadius();
        endcapRmax = endcapBounds->outerRadius();
        endcapZmin
            = fabs(endcapVolume->center().z()) - endcapBounds->halflengthZ();
        endcapZmax
            = fabs(endcapVolume->center().z()) + endcapBounds->halflengthZ();
        ACTS_VERBOSE("Outer Endcap bounds provided by configuration, "
                     "rMin/rMax/zMin/zMax = "
                     << endcapRmin
                     << ", "
                     << endcapRmax
                     << ", "
                     << endcapZmin
                     << ", "
                     << endcapZmax);
      }
      // now set the wrapping condition
      // wrapping condition can only be set if there's an inside volume
      if (insideVolume) {
        if (endcapVolume && endcapZmin < insideVolumeZmax)
          wrappingCondition = 1;
        else
          wrappingCondition = 2;
      }
    } else {
      // if no containing volumes are provided calculate the bounds from the
      // layer configuration
      // wrapping condition can only be set if there's an inside volume
      if (insideVolume) {
        if (insideVolumeRmax > volumeRmin) {
          // we need to bail out, the given volume does not fit around the other
          ACTS_ERROR("Given layer dimensions do not fit around the provided "
                     "inside volume. Bailing out."
                     << "insideVolumeRmax: "
                     << insideVolumeRmax
                     << " layerRmin: "
                     << layerRmin);
          // cleanup the memory
          negativeLayers.clear();
          centralLayers.clear();
          positiveLayers.clear();
          // return a null pointer, upstream builder will have to understand
          // this
          return nullptr;
        }
        // obvious settings
        // the new barrel inner radius is the inside volume outer radius
        barrelRmin = insideVolumeRmax;
        // the new barrel outer radius is the estimated maximal outer radius
        barrelRmax = volumeRmax;
        // regardless if the endcap setup exists or not
        endcapRmax = volumeRmax;
        endcapZmax = volumeZmax;
        // if the endcap layers exist and are within the inside volume extend
        // this is radial wrapping
        if (pLayerSetup && pLayerSetup.zBoundaries.first < insideVolumeZmax) {
          // set the barrel / endcap z division in the middle
          barrelZmax = 0.5 * (cLayerSetup.zBoundaries.second
                              + pLayerSetup.zBoundaries.first);
          // set the endcap inner radius to wrap the inner volume
          endcapRmin = insideVolumeRmax;
          // set the wrapping condition
          wrappingCondition = 1;
        } else {
          // set the barrel parameters first to the volume Zmax
          barrelZmax = volumeZmax;
          // adapt in case endcaps eixt
          if (pLayerSetup) {
            // the barrel z extend is either set to the inside z extend
            /// or into the middle of the two
            barrelZmax = cLayerSetup.zBoundaries.second < insideVolumeZmax
                ? insideVolumeZmax
                : 0.5 * (cLayerSetup.zBoundaries.second
                         + pLayerSetup.zBoundaries.first);
            // set the endcap parameters
            endcapRmin = insideVolumeRmin;
          }
          // set the wrapping condition
          wrappingCondition = 2;
        }
        // consequent setting (regardless if endcaps exist or not )
        endcapZmin = barrelZmax;

      } else {
        // no inside volume is given, wrapping conditions remains 0
        barrelRmin = volumeRmin;
        barrelRmax = volumeRmax;
        barrelZmax = pLayerSetup
            ? 0.5 * (cLayerSetup.zBoundaries.second
                     + pLayerSetup.zBoundaries.first)
            : volumeZmax;
        // endcap parameters
        endcapRmin = volumeRmin;
        endcapRmax = volumeRmax;
        endcapZmin = barrelZmax;
        endcapZmax = volumeZmax;
      }
    }  // else - no volume bounds given from translation

    // the barrel is created
    barrel = m_cfg.trackingVolumeHelper->createTrackingVolume(
        centralLayers,
        m_cfg.volumeMaterial,
        barrelRmin,
        barrelRmax,
        -barrelZmax,
        barrelZmax,
        m_cfg.volumeName + "::Barrel");

    // the negative endcap is created
    nEndcap = !negativeLayers.empty()
        ? m_cfg.trackingVolumeHelper->createTrackingVolume(
              negativeLayers,
              m_cfg.volumeMaterial,
              endcapRmin,
              endcapRmax,
              -endcapZmax,
              -endcapZmin,
              m_cfg.volumeName + "::NegativeEndcap")
        : nullptr;

    // the positive endcap is created
    pEndcap = !positiveLayers.empty()
        ? m_cfg.trackingVolumeHelper->createTrackingVolume(
              positiveLayers,
              m_cfg.volumeMaterial,
              endcapRmin,
              endcapRmax,
              endcapZmin,
              endcapZmax,
              m_cfg.volumeName + "::PositiveEndcap")
        : nullptr;

    // no wrapping condition
    if (wrappingCondition == 0) {
      // we have endcap volumes
      if (nEndcap && pEndcap) {
        // a new barrel sector
        volume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
            {nEndcap, barrel, pEndcap});
      } else  // just take the barrel as the return value
        volume = barrel;

    } else if (wrappingCondition == 1) {
      // a new barrel sector
      volume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
          {nEndcap, barrel, pEndcap});
      // now check if we need gaps as in 1
      if (fabs(insideVolumeZmax - volumeZmax) > 10e-5) {
        // create the gap volumes
        // - negative side
        nEndcap = m_cfg.trackingVolumeHelper->createGapTrackingVolume(
            m_cfg.volumeMaterial,
            insideVolumeRmax,
            volumeRmax,
            -volumeZmax,
            -barrelZmax,
            1,
            false,
            m_cfg.volumeName + "::NegativeGap");
        // - positive side
        pEndcap = m_cfg.trackingVolumeHelper->createGapTrackingVolume(
            m_cfg.volumeMaterial,
            insideVolumeRmax,
            volumeRmax,
            barrelZmax,
            volumeZmax,
            1,
            false,
            m_cfg.volumeName + "::PositiveGap");
        // update the volume with the two sides
        insideVolume
            = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
                {nEndcap, insideVolume, pEndcap});
      }
      // update the volume
      volume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
          {insideVolume, volume});

    } else if (wrappingCondition == 2) {
      // create gap volumes if needed
      if (barrelZmax > insideVolumeZmax) {
        // create the gap volumes
        auto niGap = m_cfg.trackingVolumeHelper->createGapTrackingVolume(
            m_cfg.volumeMaterial,
            insideVolumeRmin,
            volumeRmin,
            -barrelZmax,
            -insideVolumeZmax,
            1,
            false,
            m_cfg.volumeName + "::InnerNegativeGap");

        auto piGap = m_cfg.trackingVolumeHelper->createGapTrackingVolume(
            m_cfg.volumeMaterial,
            insideVolumeRmin,
            volumeRmin,
            insideVolumeZmax,
            barrelZmax,
            1,
            false,
            m_cfg.volumeName + "::InnerPositiveGap");
        // pack into a new insideVolume
        insideVolume
            = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
                {niGap, insideVolume, piGap});
      }
      // create the container of the detector
      insideVolume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
          {insideVolume, barrel});
      volume = (nEndcap && pEndcap)
          ? m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
                {nEndcap, insideVolume, pEndcap})
          : insideVolume;
    }
  } else if (outsideBounds) {
    // screen output
    ACTS_DEBUG("Building Volume without layer configuration.");
    if (insideVolume && outsideBounds) {
      // the barrel is created
      barrel = m_cfg.trackingVolumeHelper->createTrackingVolume(
          {},
          m_cfg.volumeMaterial,
          insideVolumeRmin,
          volumeRmax,
          -insideVolumeZmax,
          insideVolumeZmax,
          m_cfg.volumeName + "::Barrel");
      // pack into the appropriate container
      volume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
          {insideVolume, barrel});
      // check if necap gaps are needed
      if (fabs(insideVolumeZmax - volumeZmax) > 10e-5) {
        // the negative endcap is created
        nEndcap = m_cfg.trackingVolumeHelper->createTrackingVolume(
            {},
            m_cfg.volumeMaterial,
            insideVolumeRmin,
            volumeRmax,
            -volumeZmax,
            -insideVolumeZmax,
            m_cfg.volumeName + "::NegativeEndcap");
        // the positive endcap is created
        pEndcap = m_cfg.trackingVolumeHelper->createTrackingVolume(
            {},
            m_cfg.volumeMaterial,
            insideVolumeRmin,
            volumeRmax,
            insideVolumeZmax,
            volumeZmax,
            m_cfg.volumeName + "::PositiveEndcap");
        // pack into a the container
        volume = m_cfg.trackingVolumeHelper->createContainerTrackingVolume(
            {nEndcap, barrel, pEndcap});
      }

    } else
      volume = TrackingVolume::create(
          nullptr, outsideBounds, m_cfg.volumeMaterial);

  } else {
    ACTS_ERROR(
        "Neither layer configuration nor volume bounds given. Bailing out.");
  }
  // sign the volume
  volume->sign(GeometrySignature(m_cfg.volumeSignature));
  // now return what you have
  return volume;
}

Acts::LayerSetup
Acts::CylinderVolumeBuilder::analyzeLayerSetup(const LayerVector lVector) const
{
  // return object
  LayerSetup lSetup;
  // only if the vector is present it can actually be analyzed
  if (!lVector.empty()) {
    // we have layers
    lSetup.present = true;
    for (auto& layer : lVector) {
      // the thickness of the layer needs to be taken into account
      double thickness = layer->thickness();
      // get the center of the layer
      const Vector3D& center = layer->surfaceRepresentation().center();
      // check if it is a cylinder layer
      const CylinderLayer* cLayer
          = dynamic_cast<const CylinderLayer*>(layer.get());
      if (cLayer) {
        // set the binning to radial binning
        lSetup.binningValue = binR;
        // now we have access to all the information
        double rMinC
            = cLayer->surfaceRepresentation().bounds().r() - 0.5 * thickness;
        double rMaxC
            = cLayer->surfaceRepresentation().bounds().r() + 0.5 * thickness;
        double hZ = cLayer->surfaceRepresentation().bounds().halflengthZ();
        takeSmaller(lSetup.rBoundaries.first, rMinC);
        takeBigger(lSetup.rBoundaries.second, rMaxC);
        takeSmaller(lSetup.zBoundaries.first, center.z() - hZ);
        takeBigger(lSetup.zBoundaries.second, center.z() + hZ);
      }
      // proceed further if it is a Disc layer
      const RadialBounds* dBounds = dynamic_cast<const RadialBounds*>(
          &(layer->surfaceRepresentation().bounds()));
      if (dBounds) {
        // set the binning to radial binning
        lSetup.binningValue = binZ;
        // now we have access to all the information
        double rMinD = dBounds->rMin();
        double rMaxD = dBounds->rMax();
        double zMinD = center.z() - 0.5 * thickness;
        double zMaxD = center.z() + 0.5 * thickness;
        takeSmaller(lSetup.rBoundaries.first, rMinD);
        takeBigger(lSetup.rBoundaries.second, rMaxD);
        takeSmaller(lSetup.zBoundaries.first, zMinD);
        takeBigger(lSetup.zBoundaries.second, zMaxD);
        //!< @TODO check for Ring setup
      }
    }
  }
  return lSetup;
}
