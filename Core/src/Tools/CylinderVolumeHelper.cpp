// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderVolumeHelper.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/CylinderVolumeHelper.hpp"
#include "ACTS/Detector/GlueVolumesDescriptor.hpp"
#include "ACTS/Detector/TrackingVolume.hpp"
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Material/BinnedSurfaceMaterial.hpp"
#include "ACTS/Material/Material.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Tools/ILayerArrayCreator.hpp"
#include "ACTS/Tools/ITrackingVolumeArrayCreator.hpp"
#include "ACTS/Utilities/BinUtility.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Volumes/AbstractVolume.hpp"
#include "ACTS/Volumes/BoundarySurfaceT.hpp"
#include "ACTS/Volumes/CylinderVolumeBounds.hpp"

Acts::CylinderVolumeHelper::CylinderVolumeHelper(
    const Acts::CylinderVolumeHelper::Config& cvhConfig,
    std::unique_ptr<Logger>                   logger)
  : Acts::ITrackingVolumeHelper(), m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(cvhConfig);
}

// configuration
void
Acts::CylinderVolumeHelper::setConfiguration(
    const Acts::CylinderVolumeHelper::Config& cvhConfig)
{
  // @TODO check consistency
  // copy the configuration
  m_cfg = cvhConfig;
}

void
Acts::CylinderVolumeHelper::setLogger(std::unique_ptr<Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeHelper::createTrackingVolume(
    const LayerVector&                  layers,
    std::shared_ptr<Material>           matprop,
    std::shared_ptr<const VolumeBounds> volBounds,
    std::shared_ptr<Transform3D>        transform,
    const std::string&                  volumeName,
    BinningType                         bType) const
{
  // the final one to build / sensitive Volume / Bounds
  TrackingVolumePtr tVolume = nullptr;
  // the layer array
  std::unique_ptr<const LayerArray> layerArray = nullptr;

  // cases are:
  // (1) volBounds && transform   : use both information
  // (2) volBounds && !transform  : centered around 0, but with given bounds
  // (3) !volBounds && transform  : estimate size from layers, use transform
  // (4) !volBounds && !transform : estimate size & translation from layers
  const CylinderVolumeBounds* cylinderBounds = nullptr;
  // this is the implementation of CylinderVolumeHelper
  if (volBounds) {
    cylinderBounds = dynamic_cast<const CylinderVolumeBounds*>(volBounds.get());
    if (!cylinderBounds) {
      ACTS_WARNING(
          "[!] Problem: given bounds are not cylindrical - return nullptr");
      return tVolume;
    }
  }
  // this is only needed if layers are provided
  if (layers.size()) {
    // the raw data
    double rMinRaw = 0.;
    double rMaxRaw = 0.;
    double zMinRaw = 0.;
    double zMaxRaw = 0.;

    BinningValue bValue = binR;

    // check the dimension and fill raw data
    if (not estimateAndCheckDimension(layers,
                                      cylinderBounds,
                                      transform,
                                      rMinRaw,
                                      rMaxRaw,
                                      zMinRaw,
                                      zMaxRaw,
                                      bValue,
                                      bType)) {
      ACTS_WARNING("[!] Problem with given dimensions - return nullptr and "
                   "delete provided objects");
      // delete if newly created bounds
      if (volBounds == nullptr) delete cylinderBounds;
      return tVolume;
    }
    // get the zMin/Max
    double zMin = (transform ? transform->translation().z() : 0.)
        + (cylinderBounds ? -cylinderBounds->halflengthZ() : 0.);
    double zMax = (transform ? transform->translation().z() : 0.)
        + (cylinderBounds ? cylinderBounds->halflengthZ() : 0.);
    // get the rMin/rmAx
    double rMin = cylinderBounds ? cylinderBounds->innerRadius() : rMinRaw;
    double rMax = cylinderBounds ? cylinderBounds->outerRadius() : rMaxRaw;

    ACTS_VERBOSE("Filling the layers into an appropriate layer array - with "
                 "binningValue = "
                 << bValue);

    // create the Layer Array
    layerArray = (bValue == binR)
        ? m_cfg.layerArrayCreator->layerArray(layers, rMin, rMax, bType, bValue)
        : m_cfg.layerArrayCreator->layerArray(
              layers, zMin, zMax, bType, bValue);

  }  // layers are created and done
  // make sure the ownership of the bounds is correct
  std::shared_ptr<const VolumeBounds> volBoundsFinal = volBounds.get()
      ? volBounds
      : std::shared_ptr<const VolumeBounds>(cylinderBounds);
  // finally create the TrackingVolume
  tVolume = TrackingVolume::create(transform,
                                   volBoundsFinal,
                                   matprop,
                                   std::move(layerArray),
                                   {},
                                   {},
                                   {},
                                   volumeName);
  // screen output
  ACTS_VERBOSE(
      "Created cylindrical volume at z-position :" << tVolume->center().z());
  ACTS_VERBOSE("   created bounds : " << tVolume->volumeBounds());
  // return the constructed TrackingVolume
  return tVolume;
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeHelper::createTrackingVolume(
    const LayerVector&        layers,
    std::shared_ptr<Material> matprop,
    double                    rMin,
    double                    rMax,
    double                    zMin,
    double                    zMax,
    const std::string&        volumeName,
    BinningType               bType) const
{
  // that's what is needed
  CylinderVolumeBounds* cBounds = nullptr;

  // screen output
  ACTS_VERBOSE("Create cylindrical TrackingVolume '" << volumeName << "'.");
  ACTS_VERBOSE(
      "    -> with given dimensions of (rMin/rMax/zMin/Max) = " << rMin << " / "
                                                                << rMax
                                                                << " / "
                                                                << zMin
                                                                << " / "
                                                                << zMax);

  // check for consistency
  if (zMin > zMax || rMin > rMax) {
    ACTS_WARNING("Inconsistent dimensions given :"
                 << ((zMin > zMax) ? " zMin > zMax (" : " rMin > rMax (")
                 << ((zMin > zMax) ? zMin : rMin)
                 << " > "
                 << ((zMin > zMax) ? zMax : rMax)
                 << " ) - return 0");
    return nullptr;
  }

  // create a Transform3D and VolumeBounds out of the zMin/zMax
  double halflengthZ = 0.5 * (zMax - zMin);
  double zPosition   = 0.5 * (zMin + zMax);
  zPosition          = fabs(zPosition) < 0.1 ? 0. : zPosition;

  // now create the cylinder volume bounds
  cBounds = rMin > 0.1 ? new CylinderVolumeBounds(rMin, rMax, halflengthZ)
                       : new CylinderVolumeBounds(rMax, halflengthZ);
  // transform
  std::shared_ptr<Transform3D> transform
      = (zPosition != 0) ? std::make_shared<Transform3D>() : 0;
  if (transform) (*transform) = Translation3D(0., 0., zPosition);
  // call to the creation method with Bounds & Translation3D
  return createTrackingVolume(
      layers, matprop, VolumeBoundsPtr(cBounds), transform, volumeName, bType);
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeHelper::createGapTrackingVolume(
    std::shared_ptr<Material> matprop,
    double                    rMin,
    double                    rMax,
    double                    zMin,
    double                    zMax,
    unsigned int              materialLayers,
    bool                      cylinder,
    const std::string&        volumeName) const
{
  // screen output
  ACTS_VERBOSE("Create cylindrical gap TrackingVolume '"
               << volumeName
               << "' with (rMin/rMax/zMin/Max) = ");
  ACTS_VERBOSE('\t' << rMin << " / " << rMax << " / " << zMin << " / " << zMax);

  // assing min/max
  double min = cylinder ? rMin : zMin;
  double max = cylinder ? rMax : zMax;

  // create the layer r/z positions
  std::vector<double> layerPositions;
  if (materialLayers > 1) {
    double step = cylinder ? (max - min) / (materialLayers - 1)
                           : (max - min) / (materialLayers - 1);
    for (unsigned int il = 0; il < materialLayers; ++il)
      layerPositions.push_back(min + il * step);
  } else
    layerPositions.push_back(0.5 * (min + max));

  // now call the main method
  return createGapTrackingVolume(matprop,
                                 rMin,
                                 rMax,
                                 zMin,
                                 zMax,
                                 layerPositions,
                                 cylinder,
                                 volumeName,
                                 arbitrary);
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeHelper::createGapTrackingVolume(
    std::shared_ptr<Material>  matprop,
    double                     rMin,
    double                     rMax,
    double                     zMin,
    double                     zMax,
    const std::vector<double>& layerPositions,
    bool                       cylinder,
    const std::string&         volumeName,
    BinningType                bType) const
{
  // screen output
  ACTS_VERBOSE("Create cylindrical gap TrackingVolume '"
               << volumeName
               << "' with (rMin/rMax/zMin/Max) = ");
  ACTS_VERBOSE('\t' << rMin << " / " << rMax << " / " << zMin << " / " << zMax);

  // create the layers
  LayerVector layers;
  layers.reserve(layerPositions.size());

  std::vector<double>::const_iterator layerPropIter = layerPositions.begin();
  std::vector<double>::const_iterator layerPropEnd  = layerPositions.end();
  for (; layerPropIter != layerPropEnd; ++layerPropIter) {
    // create cylinder layers
    if (cylinder) {
      // take envelopes into account
      double zMinLayer = zMin;
      double zMaxLayer = zMax;
      // create the layer
      layers.push_back(createCylinderLayer(0.5 * (zMinLayer + zMaxLayer),
                                           (*layerPropIter),
                                           fabs(0.5 * (zMaxLayer - zMinLayer)),
                                           m_cfg.passiveLayerThickness,
                                           m_cfg.passiveLayerPhiBins,
                                           m_cfg.passiveLayerRzBins));

    } else {
      // take the envelopes into account
      double rMinLayer = rMin;
      double rMaxLayer = rMax;
      // create the layer
      layers.push_back(createDiscLayer((*layerPropIter),
                                       rMinLayer,
                                       rMaxLayer,
                                       m_cfg.passiveLayerThickness,
                                       m_cfg.passiveLayerPhiBins,
                                       m_cfg.passiveLayerRzBins));
    }
  }
  // now call the createTrackingVolume() method
  return createTrackingVolume(
      layers, matprop, rMin, rMax, zMin, zMax, volumeName, bType);
}

std::shared_ptr<const Acts::TrackingVolume>
Acts::CylinderVolumeHelper::createContainerTrackingVolume(
    const TrackingVolumeVector& volumes) const
{
  // check if you have more than one volume
  if (volumes.size() <= (size_t)1) {
    ACTS_WARNING("None (only one) TrackingVolume given to create container "
                 "volume (min required: 2) - returning 0 ");
    return nullptr;
  }
  // screen output
  std::string volumeName = "{ ";
  ACTS_VERBOSE("[start] Creating a container volume with " << volumes.size()
                                                           << " sub volumes:");
  // volumes need to be sorted in either r or z - both increasing
  // set the iterator to the volumes, the first and the end
  auto firstVolume = volumes.begin();
  auto lastVolume  = volumes.end();

  for (size_t ivol = 0; firstVolume != lastVolume; ++firstVolume, ++ivol) {
    ACTS_VERBOSE(
        "   - volume (" << ivol << ") is : " << (*firstVolume)->volumeName());
    ACTS_VERBOSE("     at position : " << (*firstVolume)->center().x() << ", "
                                       << (*firstVolume)->center().y()
                                       << ", "
                                       << (*firstVolume)->center().z());

    ACTS_VERBOSE("     with bounds : " << (*firstVolume)->volumeBounds());
    // put the name together
    volumeName += (*firstVolume)->volumeName();
    if (ivol + 1 < volumes.size()) volumeName += " | ";
  }
  // close the volume name
  volumeName += " }";
  // reset the iterator -----
  firstVolume = volumes.begin();
  --lastVolume;  // set to the last volume

  if (firstVolume == lastVolume) {
    ACTS_WARNING("Only one TrackingVolume given to create Top level volume "
                 "(min required: 2) - returning 0 ");
    return nullptr;
  }
  // get the bounds
  const CylinderVolumeBounds* firstVolumeBounds
      = dynamic_cast<const CylinderVolumeBounds*>(
          &((*firstVolume)->volumeBounds()));
  const CylinderVolumeBounds* lastVolumeBounds
      = dynamic_cast<const CylinderVolumeBounds*>(
          &((*lastVolume)->volumeBounds()));
  // check the dynamic cast
  if (!firstVolumeBounds || !lastVolumeBounds) {
    ACTS_WARNING("VolumeBounds given are not of type: CylinderVolumeBounds "
                 "(required) - returning 0 ");
    return nullptr;
  }
  // check whether it is a r-binned case or a z-binned case
  bool rCase
      = fabs(firstVolumeBounds->innerRadius() - lastVolumeBounds->innerRadius())
      > 0.1;

  // fill these ones depending on the rCase though assignment - no parsing at
  // that stage
  double zMin     = 0.;
  double zMax     = 0.;
  double rMin     = 0.;
  double rGlueMin = 0.;
  double rMax     = 0.;
  double zSep1    = 0.;
  double zSep2    = 0.;
  if (rCase) {
    zMin     = (*firstVolume)->center().z() - firstVolumeBounds->halflengthZ();
    zMax     = (*firstVolume)->center().z() + firstVolumeBounds->halflengthZ();
    zSep1    = zMin;
    zSep2    = zMax;
    rMin     = firstVolumeBounds->innerRadius();
    rGlueMin = firstVolumeBounds->outerRadius();
    rMax     = lastVolumeBounds->outerRadius();
  } else {
    zMin  = (*firstVolume)->center().z() - firstVolumeBounds->halflengthZ();
    zMax  = (*lastVolume)->center().z() + lastVolumeBounds->halflengthZ();
    zSep1 = (*firstVolume)->center().z() + firstVolumeBounds->halflengthZ();
    zSep2 = zSep1;
    rMin  = firstVolumeBounds->innerRadius();
    rMax  = firstVolumeBounds->outerRadius();
  }
  // estimate the z - position
  double zPos = 0.5 * (zMin + zMax);
  // create the HEP transform from the stuff known so far
  std::shared_ptr<Transform3D> topVolumeTransform
      = fabs(zPos) > 0.1 ? std::make_shared<Transform3D>() : nullptr;
  if (topVolumeTransform) (*topVolumeTransform) = Translation3D(0., 0., zPos);
  // create the bounds from the information gathered so far
  CylinderVolumeBounds* topVolumeBounds = fabs(rMin) > 0.1
      ? new CylinderVolumeBounds(rMin, rMax, 0.5 * fabs(zMax - zMin))
      : new CylinderVolumeBounds(rMax, 0.5 * fabs(zMax - zMin));

  // some screen output
  ACTS_VERBOSE("Container voume bounds are " << (*topVolumeBounds));

  // create the volume array with the ITrackingVolumeArrayCreator
  std::shared_ptr<const TrackingVolumeArray> volumeArray = (rCase)
      ? m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(volumes, binR)
      : m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(volumes, binZ);
  if (volumeArray == nullptr) {
    ACTS_WARNING(
        "Creation of TrackingVolume array did not succeed - returning 0 ");
    delete topVolumeBounds;
    return nullptr;
  }
  // we have the bounds and the volume array, create the volume
  std::shared_ptr<const TrackingVolume> topVolume
      = TrackingVolume::create(topVolumeTransform,
                               VolumeBoundsPtr(topVolumeBounds),
                               volumeArray,
                               volumeName);
  // glueing section
  // --------------------------------------------------------------------------------------
  if (not interGlueTrackingVolume(
          topVolume, rCase, rMin, rGlueMin, rMax, zSep1, zSep2)) {
    ACTS_WARNING("Problem with inter-glueing of TrackingVolumes (needed) - "
                 "returning 0 ");
    return nullptr;
  }

  ACTS_VERBOSE(
      "[ end ] return newly created container : " << topVolume->volumeName());

  return topVolume;
}

/** private helper method to estimate and check the dimensions of a tracking
 * volume */
bool
Acts::CylinderVolumeHelper::estimateAndCheckDimension(
    const LayerVector&            layers,
    const CylinderVolumeBounds*&  cylinderVolumeBounds,
    std::shared_ptr<Transform3D>& transform,
    double&                       rMinClean,
    double&                       rMaxClean,
    double&                       zMinClean,
    double&                       zMaxClean,
    BinningValue&                 bValue,
    BinningType) const
{
  // some verbose output

  ACTS_VERBOSE("Parsing the " << layers.size()
                              << " layers to gather overall dimensions");
  if (cylinderVolumeBounds)
    ACTS_DEBUG("Cylinder volume bounds are given: (rmin/rmax/dz) = "
               << "("
               << cylinderVolumeBounds->innerRadius()
               << "/"
               << cylinderVolumeBounds->outerRadius()
               << "/"
               << cylinderVolumeBounds->halflengthZ()
               << ")");

  // prepare for parsing the layers
  double layerRmin = 10e10;
  double layerRmax = 0.;
  double layerZmin = 10e10;
  double layerZmax = -10e10;
  bool   radial    = false;

  rMinClean = 10e10;
  rMaxClean = 0.;
  zMinClean = 10e10;
  zMaxClean = -10e10;

  // find out what is there
  for (auto& layerIter : layers) {
    // initialize
    double currentRmin = 0.;
    double currentRmax = 0.;
    double currentZmin = 0.;
    double currentZmax = 0.;
    // dynamic cast the bounds either to CylinderBounds or DiscBounds
    const CylinderBounds* cylBounds = dynamic_cast<const CylinderBounds*>(
        &(layerIter->surfaceRepresentation()).bounds());
    // cylinder bounds
    if (cylBounds) {
      radial = true;
      // get the raw data
      double currentR = cylBounds->r();
      double centerZ  = (layerIter->surfaceRepresentation()).center().z();
      // check for min/max in the cylinder bounds case
      currentRmin = currentR - (0.5 * (layerIter)->thickness());
      currentRmax = currentR + (0.5 * (layerIter)->thickness());
      currentZmin = centerZ - cylBounds->halflengthZ();
      currentZmax = centerZ + cylBounds->halflengthZ();
    }
    // dynamic cast to the DiscBounds
    const RadialBounds* discBounds = dynamic_cast<const RadialBounds*>(
        &(layerIter->surfaceRepresentation()).bounds());
    if (discBounds) {
      // check for min/max in the cylinder bounds case
      double centerZ = (layerIter->surfaceRepresentation()).center().z();
      currentRmin    = discBounds->rMin();
      currentRmax    = discBounds->rMax();
      currentZmin    = centerZ - (0.5 * (layerIter)->thickness());
      currentZmax    = centerZ + (0.5 * (layerIter)->thickness());
    }
    // the raw data
    takeSmaller(rMinClean, currentRmin);
    takeBigger(rMaxClean, currentRmax);
    takeSmaller(zMinClean, currentZmin);
    takeBigger(zMaxClean, currentZmax);
    // assign if they overrule the minima/maxima (with layers thicknesses)
    takeSmaller(layerRmin, currentRmin);
    takeBigger(layerRmax, currentRmax);
    takeSmaller(layerZmin, currentZmin);
    takeBigger(layerZmax, currentZmax);
  }

  // set the binning value
  bValue = radial ? binR : binZ;

  ACTS_VERBOSE("Estimate/check CylinderVolumeBounds from/w.r.t. enclosed "
               "layers + envelope covers");
  // the z from the layers w and w/o envelopes
  double zEstFromLayerEnv    = 0.5 * ((layerZmax) + (layerZmin));
  double halflengthFromLayer = 0.5 * fabs((layerZmax) - (layerZmin));

  bool concentric = (zEstFromLayerEnv * zEstFromLayerEnv < 0.001);

  // no CylinderBounds and Translation given - make it
  if (!cylinderVolumeBounds && !transform) {
    // create the CylinderBounds from parsed layer inputs
    cylinderVolumeBounds
        = new CylinderVolumeBounds(layerRmin, layerRmax, halflengthFromLayer);
    // and the transform
    transform = concentric ? std::make_shared<Transform3D>() : 0;
    if (transform) (*transform) = Translation3D(0., 0., zEstFromLayerEnv);
  } else if (cylinderVolumeBounds && !transform && !concentric) {
    transform    = std::make_shared<Transform3D>();
    (*transform) = Translation3D(0., 0., zEstFromLayerEnv);
  } else if (transform && !cylinderVolumeBounds) {
    // create the CylinderBounds from parsed layer inputs
    double halflengthFromLayer = 0.5 * fabs((layerZmax) - (layerZmin));
    cylinderVolumeBounds
        = new CylinderVolumeBounds(layerRmin, layerRmax, halflengthFromLayer);
  }

  ACTS_VERBOSE("    -> dimensions from layers   (rMin/rMax/zMin/zMax) = "
               << layerRmin
               << " / "
               << layerRmax
               << " / "
               << layerZmin
               << " / "
               << layerZmax);

  double zFromTransform = transform ? transform->translation().z() : 0.;
  ACTS_VERBOSE("    -> while created bounds are (rMin/rMax/zMin/zMax) = "
               << cylinderVolumeBounds->innerRadius()
               << " / "
               << cylinderVolumeBounds->outerRadius()
               << " / "
               << zFromTransform - cylinderVolumeBounds->halflengthZ()
               << " / "
               << zFromTransform + cylinderVolumeBounds->halflengthZ());

  // both is NOW given --- check it -----------------------------
  if (cylinderVolumeBounds) {
    // only check
    if (zFromTransform - cylinderVolumeBounds->halflengthZ() <= layerZmin
        && zFromTransform + cylinderVolumeBounds->halflengthZ() >= layerZmax
        && cylinderVolumeBounds->innerRadius() <= layerRmin
        && cylinderVolumeBounds->outerRadius() >= layerRmax)
      return true;
    else {
      ACTS_WARNING(
          "Provided layers are not contained by volume ! Bailing out. ");
      ACTS_WARNING("- zFromTransform: " << zFromTransform);
      ACTS_WARNING("- volumeZmin:"
                   << zFromTransform - cylinderVolumeBounds->halflengthZ()
                   << ", layerZmin: "
                   << layerZmin);
      ACTS_WARNING("- volumeZmax: "
                   << zFromTransform + cylinderVolumeBounds->halflengthZ()
                   << ", layerZmax: "
                   << layerZmax);
      ACTS_WARNING("- volumeRmin: " << cylinderVolumeBounds->innerRadius()
                                    << ", layerRmin: "
                                    << layerRmin);
      ACTS_WARNING("- volumeRmax: " << cylinderVolumeBounds->outerRadius()
                                    << ", layerRmax: "
                                    << layerRmax);
      return false;
    }
  }

  ACTS_VERBOSE("Created/Checked " << *cylinderVolumeBounds);
  return true;
}

bool
Acts::CylinderVolumeHelper::interGlueTrackingVolume(
    std::shared_ptr<const TrackingVolume> tVolume,
    bool                                  rBinned,
    double                                rMin,
    double                                rGlueMin,
    double                                rMax,
    double                                zMin,
    double                                zMax) const
{
  ACTS_VERBOSE("Glue contained TrackingVolumes of container '"
               << tVolume->volumeName()
               << "'.");

  // only go on if you have confinedVolumes
  if (tVolume->confinedVolumes()) {
    // get the glueVolumes descriptor of the top volume to register the outside
    // volumes
    const GlueVolumesDescriptor& glueDescr = tVolume->glueVolumesDescriptor();

    // now retrieve the volumes
    auto& volumes = tVolume->confinedVolumes()->arrayObjects();

    // list the volume names:
    //  and make the screen output readable
    size_t ivol = 0;
    for (auto& vol : volumes)
      ACTS_VERBOSE("[" << ivol++ << "] - volume : " << vol->volumeName());

    // the needed iterators
    auto tVolIter  = volumes.begin();
    auto tVolFirst = volumes.begin();
    auto tVolLast  = volumes.end();
    --tVolLast;
    auto tVolEnd = volumes.end();

    // the glue volumes for the description
    TrackingVolumeVector glueVolumesInnerTube;
    TrackingVolumeVector glueVolumesOuterTube;
    TrackingVolumeVector glueVolumesNegativeFace;
    TrackingVolumeVector glueVolumesPositiveFace;
    // reset ivol counter
    ivol = 0;
    // volumes of increasing r
    if (rBinned) {
      // loop over the volumes -------------------------------
      for (; tVolIter != tVolEnd;) {
        // screen output
        ACTS_VERBOSE("r-binning: Processing volume [" << ivol++ << "]");
        // for the first one
        if (tVolIter == tVolFirst) {
          addFaceVolumes((*tVolIter), tubeInnerCover, glueVolumesInnerTube);
        }
        // add this or the subvolumes to the negativeFace and positiveFace
        addFaceVolumes((*tVolIter), negativeFaceXY, glueVolumesNegativeFace);
        addFaceVolumes((*tVolIter), positiveFaceXY, glueVolumesPositiveFace);
        if (tVolIter == tVolLast) {
          addFaceVolumes((*tVolIter), tubeOuterCover, glueVolumesOuterTube);
          ++tVolIter;
        } else {
          std::shared_ptr<const TrackingVolume> tVol1 = (*tVolIter);
          std::shared_ptr<const TrackingVolume> tVol2 = (*(++tVolIter));
          glueTrackingVolumes(tVol1,
                              tubeOuterCover,
                              tVol2,
                              tubeInnerCover,
                              rMin,
                              rGlueMin,
                              rMax,
                              zMin,
                              zMax);
        }
      }
    } else {
      // volumes in increasing z
      // loop over the volumes
      for (; tVolIter != tVolEnd;) {
        // screen output
        ACTS_VERBOSE("z-binning: Processing volume '"
                     << (*tVolIter)->volumeName()
                     << "'.");
        if (tVolIter == tVolFirst)
          addFaceVolumes((*tVolIter), negativeFaceXY, glueVolumesNegativeFace);
        addFaceVolumes((*tVolIter), tubeInnerCover, glueVolumesInnerTube);
        addFaceVolumes((*tVolIter), tubeOuterCover, glueVolumesOuterTube);
        if (tVolIter == tVolLast) {
          addFaceVolumes((*tVolIter), positiveFaceXY, glueVolumesPositiveFace);
          ++tVolIter;
        } else {
          std::shared_ptr<const TrackingVolume> tVol1 = (*tVolIter);
          std::shared_ptr<const TrackingVolume> tVol2 = (*(++tVolIter));
          glueTrackingVolumes(tVol1,
                              positiveFaceXY,
                              tVol2,
                              negativeFaceXY,
                              rMin,
                              rGlueMin,
                              rMax,
                              zMin,
                              zMax);
        }
      }
    }
    // create BinnedArraysand register then to the glue volume descriptor for
    // upstream glueing
    if (glueVolumesNegativeFace.size()) {
      // create the outside volume array
      std::shared_ptr<const TrackingVolumeArray> glueVolumesNegativeFaceArray
          = m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(
              glueVolumesNegativeFace, binR);
      // register the glue voluems
      glueDescr.registerGlueVolumes(negativeFaceXY,
                                    glueVolumesNegativeFaceArray);
    }
    if (glueVolumesPositiveFace.size()) {
      // create the outside volume array
      std::shared_ptr<const TrackingVolumeArray> glueVolumesPositiveFaceArray
          = m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(
              glueVolumesPositiveFace, binR);
      // register the glue voluems
      glueDescr.registerGlueVolumes(positiveFaceXY,
                                    glueVolumesPositiveFaceArray);
    }
    if (glueVolumesInnerTube.size()) {
      // create the outside volume array
      std::shared_ptr<const TrackingVolumeArray> glueVolumesInnerTubeArray
          = m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(
              glueVolumesInnerTube, binZ);
      // register the glue voluems
      glueDescr.registerGlueVolumes(tubeInnerCover, glueVolumesInnerTubeArray);
    }
    if (glueVolumesOuterTube.size()) {
      // create the outside volume array
      std::shared_ptr<const TrackingVolumeArray> glueVolumesOuterTubeArray
          = m_cfg.trackingVolumeArrayCreator->trackingVolumeArray(
              glueVolumesOuterTube, binZ);
      // register the glue voluems
      glueDescr.registerGlueVolumes(tubeOuterCover, glueVolumesOuterTubeArray);
    }

    ACTS_VERBOSE("[GV] Register " << glueVolumesNegativeFace.size()
                                  << " volumes at face "
                                  << negativeFaceXY
                                  << ":");
    for (tVolIter = glueVolumesNegativeFace.begin();
         tVolIter != glueVolumesNegativeFace.end();
         ++tVolIter)
      ACTS_VERBOSE("   -> volume '" << (*tVolIter)->volumeName());
    ACTS_VERBOSE("[GV] Register " << glueVolumesPositiveFace.size()
                                  << " volumes at face "
                                  << positiveFaceXY);
    for (tVolIter = glueVolumesPositiveFace.begin();
         tVolIter != glueVolumesPositiveFace.end();
         ++tVolIter)
      ACTS_VERBOSE("   -> volume '" << (*tVolIter)->volumeName());
    ACTS_VERBOSE("[GV] Register " << glueVolumesInnerTube.size()
                                  << " volumes at face "
                                  << tubeInnerCover);
    for (tVolIter = glueVolumesInnerTube.begin();
         tVolIter != glueVolumesInnerTube.end();
         ++tVolIter)
      ACTS_VERBOSE("   -> volume '" << (*tVolIter)->volumeName());
    ACTS_VERBOSE("[GV] Register " << glueVolumesOuterTube.size()
                                  << " volumes at face "
                                  << tubeOuterCover);
    for (tVolIter = glueVolumesOuterTube.begin();
         tVolIter != glueVolumesOuterTube.end();
         ++tVolIter)
      ACTS_VERBOSE("   -> volume '" << (*tVolIter)->volumeName());
  }
  // return success
  return true;
}

/** private helper method to fill the glue volumes (or the volume itself in) */
void
Acts::CylinderVolumeHelper::glueTrackingVolumes(
    std::shared_ptr<const TrackingVolume> tvolOne,
    BoundarySurfaceFace                   faceOne,
    std::shared_ptr<const TrackingVolume> tvolTwo,
    BoundarySurfaceFace                   faceTwo,
    double                                rMin,
    double                                rGlueMin,
    double                                rMax,
    double                                zMin,
    double                                zMax) const
{
  // get the two gluevolume descriptors
  const GlueVolumesDescriptor& gvDescriptorOne
      = tvolOne->glueVolumesDescriptor();
  const GlueVolumesDescriptor& gvDescriptorTwo
      = tvolTwo->glueVolumesDescriptor();

  size_t volOneGlueVols = gvDescriptorOne.glueVolumes(faceOne)
      ? gvDescriptorOne.glueVolumes(faceOne)->arrayObjects().size()
      : 0;
  ACTS_VERBOSE("GlueVolumeDescriptor of volume '" << tvolOne->volumeName()
                                                  << "' has "
                                                  << volOneGlueVols
                                                  << " @ "
                                                  << faceOne);
  size_t volTwoGlueVols = gvDescriptorTwo.glueVolumes(faceTwo)
      ? gvDescriptorTwo.glueVolumes(faceTwo)->arrayObjects().size()
      : 0;
  ACTS_VERBOSE("GlueVolumeDescriptor of volume '" << tvolTwo->volumeName()
                                                  << "' has "
                                                  << volTwoGlueVols
                                                  << " @ "
                                                  << faceTwo);

  // they could still be a container though - should not happen usually
  TrackingVolumePtr glueVolOne = volOneGlueVols
      ? gvDescriptorOne.glueVolumes(faceOne)->arrayObjects()[0]
      : tvolOne;
  TrackingVolumePtr glueVolTwo = volTwoGlueVols
      ? gvDescriptorTwo.glueVolumes(faceTwo)->arrayObjects()[0]
      : tvolTwo;
  // check the cases
  if (volOneGlueVols <= 1 && volTwoGlueVols <= 1) {
    // (i) one -> one
    ACTS_VERBOSE("      glue : one[ " << glueVolOne->volumeName() << " @ "
                                      << faceOne
                                      << " ]-to-one[ "
                                      << glueVolTwo->volumeName()
                                      << " @ "
                                      << faceTwo
                                      << " ]");
    // one to one is easy
    glueVolOne->glueTrackingVolume(faceOne, glueVolTwo, faceTwo);

  } else if (volOneGlueVols <= 1) {
    // (ii) one -> many
    ACTS_VERBOSE("      glue : one[ " << glueVolOne->volumeName() << " @ "
                                      << faceOne
                                      << " ]-to-many[ "
                                      << tvolTwo->volumeName()
                                      << " @ "
                                      << faceTwo
                                      << " ]");
    glueVolOne->glueTrackingVolumes(
        faceOne, gvDescriptorTwo.glueVolumes(faceTwo), faceTwo);
  } else if (volTwoGlueVols <= 1) {
    // (iii) many -> one
    ACTS_VERBOSE("      glue : many[ " << tvolOne->volumeName() << " @ "
                                       << faceOne
                                       << " ]-to-one[ "
                                       << glueVolTwo->volumeName()
                                       << " @ "
                                       << faceTwo
                                       << " ]");
    glueVolTwo->glueTrackingVolumes(
        faceTwo, gvDescriptorOne.glueVolumes(faceOne), faceOne);
  } else {
    // (iv) glue array to array
    ACTS_VERBOSE("      glue : many[ " << tvolOne->volumeName() << " @ "
                                       << faceOne
                                       << " ]-to-many[ "
                                       << tvolTwo->volumeName()
                                       << " @ "
                                       << faceTwo
                                       << " ]");

    // create the BoundarySurface as shared pointer
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> boundarySurface
        = nullptr;

    // the transform of the new boundary surface
    std::shared_ptr<Transform3D> transform = nullptr;
    if (fabs(zMin + zMax) > 0.1) {
      // it's not a concentric cylinder, so create a transform
      auto pTransform = std::make_shared<Transform3D>();
      (*pTransform)   = Translation3D(Vector3D(0., 0., 0.5 * (zMin + zMax)));
      transform       = pTransform;
    }
    // 2 cases: r-Binning and zBinning
    if (faceOne == cylinderCover || faceOne == tubeOuterCover) {
      // (1) create the BoundaryCylinderSurface
      // now create the CylinderSurface
      std::unique_ptr<const Surface> cSurface(
          new CylinderSurface(transform, rGlueMin, 0.5 * (zMax - zMin)));
      ACTS_VERBOSE("             creating a new cylindrical boundary surface "
                   "with bounds = "
                   << cSurface->bounds());
      boundarySurface
          = std::make_shared<const BoundarySurfaceT<TrackingVolume>>(
              std::move(cSurface),
              gvDescriptorOne.glueVolumes(faceOne),
              gvDescriptorTwo.glueVolumes(faceTwo));
    } else {
      // (2) create the BoundaryDiscSurface, in that case the zMin/zMax provided
      // are both the position of the disk in question
      std::unique_ptr<const Surface> dSurface(
          new DiscSurface(transform, rMin, rMax));
      ACTS_VERBOSE("             creating a new disc-like boundary surface "
                   "with bounds = "
                   << dSurface->bounds());
      boundarySurface
          = std::make_shared<const BoundarySurfaceT<TrackingVolume>>(
              std::move(dSurface),
              gvDescriptorOne.glueVolumes(faceOne),
              gvDescriptorTwo.glueVolumes(faceTwo));
    }
    // update the volume with the boundary surface accordingly
    // it's safe to access directly, they can not be nullptr
    for (auto& oneVolume : gvDescriptorOne.glueVolumes(faceOne)->arrayObjects())
      oneVolume->updateBoundarySurface(faceOne, boundarySurface);
    for (auto& twoVolume : gvDescriptorTwo.glueVolumes(faceTwo)->arrayObjects())
      twoVolume->updateBoundarySurface(faceTwo, boundarySurface);
  }  // end of case (iv)
}

/** Private method - helper method not to duplicate code */
void
Acts::CylinderVolumeHelper::addFaceVolumes(
    std::shared_ptr<const TrackingVolume> tvol,
    BoundarySurfaceFace                   glueFace,
    TrackingVolumeVector&                 vols) const
{
  ACTS_VERBOSE("Adding face volumes of face " << glueFace << " for the volume '"
                                              << tvol->volumeName()
                                              << "'.");
  // retrieve the gluevolume descriptor
  const GlueVolumesDescriptor& gvDescriptor = tvol->glueVolumesDescriptor();
  // if volumes are registered: take them
  if (gvDescriptor.glueVolumes(glueFace)) {
    // get the navigation level subvolumes
    auto volIter = gvDescriptor.glueVolumes(glueFace)->arrayObjects().begin();
    auto volEnd  = gvDescriptor.glueVolumes(glueFace)->arrayObjects().end();
    for (; volIter != volEnd; ++volIter) {
      ACTS_VERBOSE("   -> adding : " << (*volIter)->volumeName());
      vols.push_back(*volIter);
    }
    // screen output
    ACTS_VERBOSE(vols.size()
                 << " navigation volumes registered as glue volumes.");
  } else {
    // the volume itself is on navigation level
    ACTS_VERBOSE("     -> adding only volume itself (at navigation level).");
    vols.push_back(tvol);
  }
}

std::shared_ptr<const Acts::Layer>
Acts::CylinderVolumeHelper::createCylinderLayer(double z,
                                                double r,
                                                double halflengthZ,
                                                double thickness,
                                                int    binsPhi,
                                                int    binsZ) const
{
  ACTS_VERBOSE("Creating a CylinderLayer at position " << z << " and radius "
                                                       << r);
  // positioning
  std::shared_ptr<Transform3D> transform = 0;
  transform = (fabs(z) > 0.1) ? std::make_shared<Transform3D>() : 0;
  if (transform) (*transform) = Translation3D(0., 0., z);

  // z-binning
  BinUtility layerBinUtility(
      binsZ, z - halflengthZ, z + halflengthZ, open, binZ);
  if (binsPhi == 1) {
    // the BinUtility for the material
    // ---------------------> create material for the layer surface
    ACTS_VERBOSE(" -> Preparing the binned material with " << binsZ
                                                           << " bins in Z. ");

  } else {  // break the phi symmetry
    // update the BinUtility: local position on Cylinder is rPhi, z
    BinUtility layerBinUtilityPhiZ(
        binsPhi, -r * M_PI, +r * M_PI, closed, binPhi);
    layerBinUtilityPhiZ += layerBinUtility;
    // ---------------------> create material for the layer surface
    ACTS_VERBOSE(
        " -> Preparing the binned material with " << binsPhi << " / " << binsZ
                                                  << " bins in phi / Z. ");
  }
  // @TODO create the SurfaceMaterial
  // bounds for cylinderical surface
  CylinderBounds* cylinderBounds = new CylinderBounds(r, halflengthZ);
  // create the cylinder
  return CylinderLayer::create(
      transform,
      std::shared_ptr<const CylinderBounds>(cylinderBounds),
      nullptr,
      thickness);
}

std::shared_ptr<const Acts::Layer>
Acts::CylinderVolumeHelper::createDiscLayer(double z,
                                            double rMin,
                                            double rMax,
                                            double thickness,
                                            int    binsPhi,
                                            int    binsR) const
{
  ACTS_VERBOSE("Creating a DiscLayer at position " << z << " and rMin/rMax "
                                                   << rMin
                                                   << " / "
                                                   << rMax);

  // positioning
  std::shared_ptr<Transform3D> transform
      = fabs(z) > 0.1 ? std::make_shared<Transform3D>() : 0;
  if (transform) (*transform) = Translation3D(0., 0., z);

  // R is the primary binning for the material
  BinUtility materialBinUtility(binsR, rMin, rMax, open, binR);
  if (binsPhi == 1) {
    ACTS_VERBOSE(" -> Preparing the binned material with " << binsR
                                                           << " bins in R. ");
  } else {
    // also binning in phi chosen
    materialBinUtility += BinUtility(binsPhi, -M_PI, M_PI, closed, binPhi);
    ACTS_VERBOSE(
        " -> Preparing the binned material with " << binsPhi << " / " << binsR
                                                  << " bins in phi / R. ");
  }

  // @TODO create the SurfaceMaterial
  // bounds for disk-like surface
  RadialBounds* discBounds = new RadialBounds(rMin, rMax);
  // create the disc
  return DiscLayer::create(transform,
                           std::shared_ptr<const DiscBounds>(discBounds),
                           nullptr,
                           thickness);
}
