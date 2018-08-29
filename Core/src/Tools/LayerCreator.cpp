// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// LayerCreator.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Tools/LayerCreator.hpp"
#include <cmath>
#include <set>
#include "Acts/Layers/CylinderLayer.hpp"
#include "Acts/Layers/DiscLayer.hpp"
#include "Acts/Layers/ProtoLayer.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

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
Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>& surfaces,
                                  size_t                             binsPhi,
                                  size_t                             binsZ,
                                  boost::optional<ProtoLayer> _protoLayer,
                                  std::shared_ptr<const Transform3D>  transform,
                                  std::unique_ptr<ApproachDescriptor> ad) const
{

  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  // remaining layer parameters
  double layerR
      = 0.5 * (protoLayer.minR - protoLayer.envR.first + protoLayer.maxR
               + protoLayer.envR.second);
  double binPosZ = 0.5 * (protoLayer.minZ + protoLayer.maxZ);

  double envZShift = 0.5 * (-protoLayer.envZ.first + protoLayer.envZ.second);
  double layerZ    = binPosZ + envZShift;
  double layerHalfZ
      = std::abs(protoLayer.maxZ + protoLayer.envZ.second - layerZ);
  double layerThickness = (protoLayer.maxR - protoLayer.minR)
      + protoLayer.envR.first + protoLayer.envR.second;

  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R     = " << layerR);
  ACTS_VERBOSE(" - from R min/max   = " << protoLayer.minR << " / "
                                        << protoLayer.maxR);
  ACTS_VERBOSE(" - with R thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envR.first << " / "
                                        << protoLayer.envR.second);

  ACTS_VERBOSE(" - with z min/max   = " << protoLayer.minZ << " (-"
                                        << protoLayer.envZ.first
                                        << ") / "
                                        << protoLayer.maxZ
                                        << " (+"
                                        << protoLayer.envZ.second
                                        << ")");

  ACTS_VERBOSE(" - z center         = " << layerZ);

  // create the layer transforms if not given
  // we need to transform in case layerZ != 0, so that the layer will be
  // correctly defined using the halflength
  if (!transform) {
    // double shift = -(layerZ + envZShift);
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., -layerZ));
    ACTS_VERBOSE(" - layer z shift  = " << -layerZ);
  }

  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.minPhi << " / "
                                        << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size() << " ordered in ( "
                                        << binsPhi
                                        << " x "
                                        << binsZ
                                        << ")");
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
        surfaces, binsPhi, binsZ, protoLayer, nullptr);

    checkBinning(*sArray);
  }

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
Acts::LayerCreator::cylinderLayer(const std::vector<const Surface*>& surfaces,
                                  BinningType                        bTypePhi,
                                  BinningType                        bTypeZ,
                                  boost::optional<ProtoLayer> _protoLayer,
                                  std::shared_ptr<const Transform3D>  transform,
                                  std::unique_ptr<ApproachDescriptor> ad) const
{

  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  // remaining layer parameters
  double layerR
      = 0.5 * (protoLayer.minR - protoLayer.envR.first + protoLayer.maxR
               + protoLayer.envR.second);
  double binPosZ   = 0.5 * (protoLayer.minZ + protoLayer.maxZ);
  double envZShift = 0.5 * (-protoLayer.envZ.first + protoLayer.envZ.second);
  double layerZ    = binPosZ + envZShift;

  double layerHalfZ
      = std::abs(protoLayer.maxZ + protoLayer.envZ.second - layerZ);

  double layerThickness = (protoLayer.maxR - protoLayer.minR)
      + protoLayer.envR.first + protoLayer.envR.second;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a cylindrical Layer:");
  ACTS_VERBOSE(" - with layer R     = " << layerR);
  ACTS_VERBOSE(" - from R min/max   = " << protoLayer.minR << " / "
                                        << protoLayer.maxR);
  ACTS_VERBOSE(" - with R thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envR.first << " / "
                                        << protoLayer.envR.second);
  ACTS_VERBOSE(" - with z min/max   = " << protoLayer.minZ << " (-"
                                        << protoLayer.envZ.first
                                        << ") / "
                                        << protoLayer.maxZ
                                        << " (+"
                                        << protoLayer.envZ.second
                                        << ")");
  ACTS_VERBOSE(" - z center         = " << layerZ);

  // create the layer transforms if not given
  // we need to transform in case layerZ != 0, so that the layer will be
  // correctly defined using the halflength
  if (!transform && bTypeZ == equidistant) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., -layerZ));
    ACTS_VERBOSE(" - layer z shift    = " << -layerZ);
  }

  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.minPhi << " / "
                                        << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size() << "");

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnCylinder(
        surfaces, bTypePhi, bTypeZ, protoLayer, nullptr);

    checkBinning(*sArray);
  }

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
Acts::LayerCreator::discLayer(const std::vector<const Surface*>&  surfaces,
                              size_t                              binsR,
                              size_t                              binsPhi,
                              boost::optional<ProtoLayer>         _protoLayer,
                              std::shared_ptr<const Transform3D>  transform,
                              std::unique_ptr<ApproachDescriptor> ad) const
{
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  double layerZ
      = 0.5 * (protoLayer.minZ - protoLayer.envZ.first + protoLayer.maxZ
               + protoLayer.envZ.second);
  double layerThickness = (protoLayer.maxZ - protoLayer.minZ)
      + protoLayer.envZ.first + protoLayer.envZ.second;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position    = " << layerZ);
  ACTS_VERBOSE(" - from Z min/max   = " << protoLayer.minZ << " / "
                                        << protoLayer.maxZ);
  ACTS_VERBOSE(" - with Z thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envZ.first << " / "
                                        << protoLayer.envZ.second);
  ACTS_VERBOSE(" - with R min/max   = " << protoLayer.minR << " (-"
                                        << protoLayer.envR.first
                                        << ") / "
                                        << protoLayer.maxR
                                        << " (+"
                                        << protoLayer.envR.second
                                        << ")");
  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.minPhi << " / "
                                        << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules    = " << surfaces.size() << " ordered in ( "
                                       << binsR
                                       << " x "
                                       << binsPhi
                                       << ")");

  // create the layer transforms if not given
  if (!transform) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., layerZ));
  }
  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
        surfaces, binsR, binsPhi, protoLayer, transform);

    checkBinning(*sArray);
  }

  // create the share disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(
      protoLayer.minR - protoLayer.envR.first,
      protoLayer.maxR + protoLayer.envR.second);

  // create the layers
  // we use the same transform here as for the layer itself
  // for disk this is fine since we don't bin in Z, so does not matter
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
                              BinningType                         bTypeR,
                              BinningType                         bTypePhi,
                              boost::optional<ProtoLayer>         _protoLayer,
                              std::shared_ptr<const Transform3D>  transform,
                              std::unique_ptr<ApproachDescriptor> ad) const
{
  ProtoLayer protoLayer = _protoLayer ? *_protoLayer : ProtoLayer(surfaces);

  double layerZ
      = 0.5 * (protoLayer.minZ - protoLayer.envZ.first + protoLayer.maxZ
               + protoLayer.envZ.second);
  double layerThickness = std::abs(protoLayer.maxZ - protoLayer.minZ)
      + protoLayer.envZ.first + protoLayer.envZ.second;

  // adjust the layer radius
  ACTS_VERBOSE("Creating a disk Layer:");
  ACTS_VERBOSE(" - at Z position    = " << layerZ);
  ACTS_VERBOSE(" - from Z min/max   = " << protoLayer.minZ << " / "
                                        << protoLayer.maxZ);
  ACTS_VERBOSE(" - with Z thickness = " << layerThickness);
  ACTS_VERBOSE("   - incl envelope  = " << protoLayer.envZ.first << " / "
                                        << protoLayer.envZ.second);
  ACTS_VERBOSE(" - with R min/max   = " << protoLayer.minR << " (-"
                                        << protoLayer.envR.first
                                        << ") / "
                                        << protoLayer.maxR
                                        << " (+"
                                        << protoLayer.envR.second
                                        << ")");
  ACTS_VERBOSE(" - with phi min/max = " << protoLayer.minPhi << " / "
                                        << protoLayer.maxPhi);
  ACTS_VERBOSE(" - # of modules     = " << surfaces.size());

  // create the layer transforms if not given
  if (!transform) {
    transform
        = std::make_shared<const Transform3D>(Translation3D(0., 0., layerZ));
  }

  // create the surface array
  std::unique_ptr<SurfaceArray> sArray;
  if (!surfaces.empty()) {
    sArray = m_cfg.surfaceArrayCreator->surfaceArrayOnDisc(
        surfaces, bTypeR, bTypePhi, protoLayer, transform);

    checkBinning(*sArray);
  }

  // create the shared disc bounds
  auto dBounds = std::make_shared<const RadialBounds>(
      protoLayer.minR - protoLayer.envR.first,
      protoLayer.maxR + protoLayer.envR.second);

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
Acts::LayerCreator::planeLayer(
    const std::vector<const Surface*>& /**surfaces*/,
    double /**envelopeXY*/,
    double /**envelopeZ*/,
    size_t /**binsX*/,
    size_t /**binsY*/,
    const std::shared_ptr<const Transform3D>& /**transform*/,
    std::unique_ptr<ApproachDescriptor> /**ad*/) const
{
  //@todo implement
  return nullptr;
}

void
Acts::LayerCreator::associateSurfacesToLayer(Layer& layer) const
{
  if (layer.surfaceArray() != nullptr) {
    auto surfaces = layer.surfaceArray()->surfaces();

    for (auto& surface : surfaces) {
      auto mutableSurface = const_cast<Surface*>(surface);
      mutableSurface->associateLayer(layer);
    }
  }
}

bool
Acts::LayerCreator::checkBinning(const SurfaceArray& sArray) const
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
    if (!sArray.isValidBin(b)) {
      continue;
    }
    for (const auto& srf : binContent) {
      accessibleSurfaces.insert(srf);
    }
    if (binContent.empty()) {
      nEmptyBins++;
    }
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

  if (!diff.empty()) {
    ACTS_ERROR(" -- Not all sensitive surfaces are accessible through binning. "
               "sensitive: "
               << sensitiveSurfaces.size()
               << "    accessible: "
               << accessibleSurfaces.size());

    // print all inaccessibles
    ACTS_ERROR(" -- Inaccessible surfaces: ");
    for (const auto& srf : diff) {
      // have to choose BinningValue here
      Vector3D ctr = srf->binningPosition(binR);
      ACTS_ERROR(" Surface(x=" << ctr.x() << ", y=" << ctr.y() << ", z="
                               << ctr.z()
                               << ", r="
                               << ctr.perp()
                               << ", phi="
                               << LA::phi(ctr)
                               << ")");
    }

  } else {
    ACTS_VERBOSE(" -- All sensitive surfaces are accessible through binning.");
  }

  return nEmptyBins == 0 && diff.empty();
}
