// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerCreator.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Tools/LayerCreator.hpp"
#include <cmath>
#include <set>
#include "ACTS/Layers/CylinderLayer.hpp"
#include "ACTS/Layers/DiscLayer.hpp"
#include "ACTS/Layers/ProtoLayer.hpp"
#include "ACTS/Surfaces/CylinderBounds.hpp"
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/RadialBounds.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

Acts::LayerCreator::LayerCreator(const Acts::LayerCreator::Config& lcConfig,
                                 std::unique_ptr<const Logger>     logger)
  : m_cfg(), m_logger(std::move(logger))
{
  setConfiguration(lcConfig);
}

void
Acts::LayerCreator::setConfiguration(const Acts::LayerCreator::Config& lcConfig)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = lcConfig;
}

void
Acts::LayerCreator::setLogger(std::unique_ptr<const Logger> newLogger)
{
  m_logger = std::move(newLogger);
}

Acts::MutableLayerPtr
Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>&  surfaces,
                                  double                              envelopeR,
                                  double                              envelopeZ,
                                  size_t                              binsPhi,
                                  size_t                              binsZ,
                                  std::shared_ptr<const Transform3D>  transform,
                                  std::unique_ptr<ApproachDescriptor> ad) const
{

  ProtoLayer protoLayer(surfaces);

  // remaining layer parameters
  double layerR         = 0.5 * (protoLayer.minR + protoLayer.maxR);
  double layerHalfZ     = protoLayer.maxZ;
  double layerThickness = (protoLayer.maxR - protoLayer.minR) + 2 * envelopeR;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R    = " << layerR);
  ACTS_VERBOSE(" - from R min/max  = " << protoLayer.minR << " / "
                                       << protoLayer.maxR);
  ACTS_VERBOSE(" - with z min/max  = " << -layerHalfZ << " / " << layerHalfZ);
  ACTS_VERBOSE(" - with thickness  = " << (protoLayer.maxR - protoLayer.minR));
  ACTS_VERBOSE("   and tolerance   = " << envelopeR);
  ACTS_VERBOSE(" - and phi min/max = " << protoLayer.minPhi << " / "
                                       << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( "
                                       << binsPhi
                                       << " x "
                                       << binsZ
                                       << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
          surfaces, binsPhi, binsZ, protoLayer, transform);

  checkBinning(*sArray);

  // create the layer and push it back
  std::shared_ptr<const CylinderBounds> cBounds(
      new CylinderBounds(layerR, layerHalfZ + envelopeZ));

  // create the layer
  MutableLayerPtr cLayer = CylinderLayer::create(transform,
                                                 cBounds,
                                                 std::move(sArray),
                                                 layerThickness,
                                                 std::move(ad),
                                                 active);

  if (!cLayer) ACTS_ERROR("Creation of cylinder layer did not succeed!");
  associateSurfacesToLayer(*cLayer);

  // now return
  return cLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>& surfaces,
                                  double                             layerRmin,
                                  double                             layerRmax,
                                  double                             layerHalfZ,
                                  BinningType                        bTypePhi,
                                  BinningType                        bTypeZ,
                                  std::shared_ptr<const Transform3D> transform,
                                  std::unique_ptr<ApproachDescriptor> ad) const
{

  // create protolayer but overwrite specified values
  ProtoLayer protoLayer(surfaces);
  protoLayer.minR = layerRmin;
  protoLayer.maxR = layerRmax;
  protoLayer.minZ = -layerHalfZ;
  protoLayer.maxZ = layerHalfZ;

  // remaining layer parameters
  double layerR         = 0.5 * (layerRmin + layerRmax);
  double layerThickness = layerRmax - layerRmin;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R    = " << layerR);
  ACTS_VERBOSE(" - from R min/max  = " << layerRmin << " / " << layerRmax);
  ACTS_VERBOSE(" - with z min/max  = " << -layerHalfZ << " / " << layerHalfZ);
  ACTS_VERBOSE(" - with thickness  = " << layerThickness);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
          surfaces, bTypePhi, bTypeZ, protoLayer, transform);

  checkBinning(*sArray);

  // create the layer and push it back
  std::shared_ptr<const CylinderBounds> cBounds(
      new CylinderBounds(layerR, layerHalfZ));

  // create the layer
  MutableLayerPtr cLayer = CylinderLayer::create(transform,
                                                 cBounds,
                                                 std::move(sArray),
                                                 layerThickness,
                                                 std::move(ad),
                                                 active);

  if (!cLayer) ACTS_ERROR("Creation of cylinder layer did not succeed!");
  associateSurfacesToLayer(*cLayer);

  // now return
  return cLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>&  surfaces,
                                  double                              envelopeR,
                                  double                              envelopeZ,
                                  BinningType                         bTypePhi,
                                  BinningType                         bTypeZ,
                                  std::shared_ptr<const Transform3D>  transform,
                                  std::unique_ptr<ApproachDescriptor> ad) const
{
  ProtoLayer protoLayer(surfaces);

  // remaining layer parameters
  double layerR         = 0.5 * (protoLayer.minR + protoLayer.maxR);
  double layerHalfZ     = protoLayer.maxZ;
  double layerThickness = (protoLayer.maxR - protoLayer.minR) + 2 * envelopeR;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R    = " << layerR);
  ACTS_VERBOSE(" - from R min/max  = " << protoLayer.minR << " / "
                                       << protoLayer.maxR);
  ACTS_VERBOSE(" - with z min/max  = " << -layerHalfZ << " / " << layerHalfZ);
  ACTS_VERBOSE(" - with thickness  = " << (protoLayer.maxR - protoLayer.minR));
  ACTS_VERBOSE("   and tolerance   = " << envelopeR);
  ACTS_VERBOSE(" - and phi min/max = " << protoLayer.minPhi << " / "
                                       << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
          surfaces, bTypePhi, bTypeZ, protoLayer, transform);

  checkBinning(*sArray);

  // create the layer and push it back
  std::shared_ptr<const CylinderBounds> cBounds(
      new CylinderBounds(layerR, layerHalfZ + envelopeZ));

  // create the layer
  MutableLayerPtr cLayer = CylinderLayer::create(transform,
                                                 cBounds,
                                                 std::move(sArray),
                                                 layerThickness,
                                                 std::move(ad),
                                                 active);

  if (!cLayer) ACTS_ERROR("Creation of cylinder layer did not succeed!");
  associateSurfacesToLayer(*cLayer);

  // now return
  return cLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::discLayer(const std::vector<const Surface*>&  surfaces,
                              double                              envelopeMinR,
                              double                              envelopeMaxR,
                              double                              envelopeZ,
                              size_t                              binsR,
                              size_t                              binsPhi,
                              std::shared_ptr<const Transform3D>  transform,
                              std::unique_ptr<ApproachDescriptor> ad) const
{

  ProtoLayer protoLayer(surfaces);

  double layerZ         = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  double layerThickness = (protoLayer.maxZ - protoLayer.minZ) + 2 * envelopeZ;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position   = " << layerZ);
  ACTS_VERBOSE(" - from R min/max  = " << protoLayer.minR << " / "
                                       << protoLayer.maxR);
  ACTS_VERBOSE(" - with thickness  = " << (protoLayer.maxZ - protoLayer.minZ));
  ACTS_VERBOSE("   and tolerance   = " << envelopeZ);
  ACTS_VERBOSE(" - and phi min/max = " << protoLayer.minPhi << " / "
                                       << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( "
                                       << binsR
                                       << " x "
                                       << binsPhi
                                       << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
          surfaces, binsR, binsPhi, protoLayer, transform);

  checkBinning(*sArray);

  // create the share disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(
      protoLayer.minR - envelopeMinR, protoLayer.maxR + envelopeMaxR);

  // create the layer transforms if not given
  if (!transform) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., layerZ));
  }

  // create the layers
  MutableLayerPtr dLayer = DiscLayer::create(transform,
                                             dBounds,
                                             std::move(sArray),
                                             layerThickness,
                                             std::move(ad),
                                             active);

  if (!dLayer) ACTS_ERROR("Creation of disc layer did not succeed!");
  associateSurfacesToLayer(*dLayer);
  // return the layer
  return dLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::discLayer(const std::vector<const Surface*>&  surfaces,
                              double                              layerZmin,
                              double                              layerZmax,
                              double                              layerRmin,
                              double                              layerRmax,
                              BinningType                         bTypeR,
                              BinningType                         bTypePhi,
                              std::shared_ptr<const Transform3D>  transform,
                              std::unique_ptr<ApproachDescriptor> ad) const
{
  // set up proto layer but overwrite with provided values
  ProtoLayer protoLayer(surfaces);
  protoLayer.minZ = layerZmin;
  protoLayer.maxZ = layerZmax;
  protoLayer.minR = layerRmin;
  protoLayer.maxR = layerRmax;

  double layerZ         = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  double layerThickness = std::abs(protoLayer.maxZ - protoLayer.minZ);

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position   = " << layerZ);
  ACTS_VERBOSE(" - from R min/max  = " << protoLayer.minR << " / "
                                       << protoLayer.maxR);
  ACTS_VERBOSE(" - with thickness  = " << layerThickness);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
          surfaces, bTypeR, bTypePhi, protoLayer, transform);

  checkBinning(*sArray);

  // create the shared disc bounds
  auto dBounds
      = std::make_shared<const RadialBounds>(protoLayer.minR, protoLayer.maxR);

  // create the layer transforms if not given
  if (!transform) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., layerZ));
  }

  // create the layers
  MutableLayerPtr dLayer = DiscLayer::create(transform,
                                             dBounds,
                                             std::move(sArray),
                                             layerThickness,
                                             std::move(ad),
                                             active);
  if (!dLayer) ACTS_ERROR("Creation of disc layer did not succeed!");
  associateSurfacesToLayer(*dLayer);
  // return the layer
  return dLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::discLayer(const std::vector<const Surface*>&  surfaces,
                              double                              envelopeMinR,
                              double                              envelopeMaxR,
                              double                              envelopeZ,
                              BinningType                         bTypeR,
                              BinningType                         bTypePhi,
                              std::shared_ptr<const Transform3D>  transform,
                              std::unique_ptr<ApproachDescriptor> ad) const
{

  ProtoLayer protoLayer(surfaces);

  // layer parametres
  double layerZ         = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  double layerThickness = (protoLayer.maxZ - protoLayer.minZ) + 2 * envelopeZ;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position   = " << layerZ);
  ACTS_VERBOSE(" - from R min/max  = " << protoLayer.minR << " / "
                                       << protoLayer.maxR);
  ACTS_VERBOSE(" - with thickness  = " << (protoLayer.maxZ - protoLayer.minZ));
  ACTS_VERBOSE("   and tolerance   = " << envelopeZ);
  ACTS_VERBOSE(" - and phi min/max = " << protoLayer.minPhi << " / "
                                       << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << ")");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray
      = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
          surfaces, bTypeR, bTypePhi, protoLayer, transform);

  checkBinning(*sArray);

  // create the share disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(
      protoLayer.minR - envelopeMinR, protoLayer.maxR + envelopeMaxR);

  // create the layer transforms if not given
  if (!transform) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., layerZ));
  }

  // create the layers
  MutableLayerPtr dLayer = DiscLayer::create(transform,
                                             dBounds,
                                             std::move(sArray),
                                             layerThickness,
                                             std::move(ad),
                                             active);
  if (!dLayer) ACTS_ERROR("Creation of disk layer did not succeed!");
  associateSurfacesToLayer(*dLayer);
  // return the layer
  return dLayer;
}

Acts::MutableLayerPtr
Acts::LayerCreator::planeLayer(
    const std::vector<const Surface*>& /**surfaces*/,
    double /**envelopeXY*/,
    double /**envelopeZ*/,
    size_t /**binsX*/,
    size_t /**binsY*/,
    std::shared_ptr<const Transform3D> /**transform*/,
    std::unique_ptr<ApproachDescriptor> /**ad*/) const
{
  //@todo implement
  return nullptr;
}

void
Acts::LayerCreator::associateSurfacesToLayer(Layer& layer) const
{
  auto surfaces = layer.surfaceArray()->surfaces();

  for (auto& surface : surfaces) {
    auto mutableSurface = const_cast<Surface*>(surface);
    mutableSurface->associateLayer(layer);
  }
}

bool
Acts::LayerCreator::checkBinning(SurfaceArray& sArray) const
{

  // do consistency check: can we access all sensitive surfaces
  // through the binning? If not, surfaces get lost and the binning does not
  // work

  ACTS_VERBOSE("Performing consistency check")

  std::vector<const Surface*> surfaces = sArray.surfaces();
  std::set<const Surface*> sensitiveSurfaces(surfaces.begin(), surfaces.end());
  std::set<const Surface*> accessibleSurfaces;
  size_t                   nEmptyBins   = 0;
  size_t                   nBinsChecked = 0;

  // iterate over all bins
  size_t size = sArray.size();
  for (size_t b = 0; b < size; ++b) {
    std::vector<const Surface*> binContent = sArray.at(b);
    // we don't check under/overflow bins
    if (!sArray.isValidBin(b)) continue;
    for (const auto& srf : binContent) {
      accessibleSurfaces.insert(srf);
    }
    if (binContent.size() == 0) nEmptyBins++;
    nBinsChecked++;
  }

  std::vector<const Acts::Surface*> diff;
  std::set_difference(sensitiveSurfaces.begin(),
                      sensitiveSurfaces.end(),
                      accessibleSurfaces.begin(),
                      accessibleSurfaces.end(),
                      std::inserter(diff, diff.begin()));

  ACTS_VERBOSE(" - Checked " << nBinsChecked << " valid bins");

  if (nEmptyBins > 0) {
    ACTS_ERROR(" -- Not all bins point to surface. " << nEmptyBins << " empty");
  } else {
    ACTS_VERBOSE(" -- All bins point to a surface");
  }

  if (diff.size() != 0) {
    ACTS_ERROR(" -- Not all sensitive surfaces are acessible through binning. "
               "sensitive: "
               << sensitiveSurfaces.size()
               << "    accessible: "
               << accessibleSurfaces.size());
  } else {
    ACTS_VERBOSE(" -- All sensitive surfaces are accessible through binning.");
  }

  return nEmptyBins == 0 && diff.size() == 0;
}
