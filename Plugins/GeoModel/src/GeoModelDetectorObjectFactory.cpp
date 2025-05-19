// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConverters.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"

#include <algorithm>
#include <iostream>
#include <typeinfo>

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoPcon.h>
#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoShapeSubtraction.h>
#include <GeoModelKernel/GeoShapeUnion.h>
#include <GeoModelKernel/GeoSimplePolygonBrep.h>
#include <GeoModelKernel/GeoTrd.h>
#include <GeoModelKernel/GeoTube.h>
#include <GeoModelKernel/GeoTubs.h>
#include <GeoModelHelpers/GeoShapeUtils.h>


namespace Acts {
GeoModelDetectorObjectFactory::GeoModelDetectorObjectFactory(
    const Config &cfg, std::unique_ptr<const Logger> mlogger)
    : m_logger(std::move(mlogger)), m_cfg(cfg) {}

void GeoModelDetectorObjectFactory::construct(
    Cache &cache, const GeometryContext &gctx, const GeoModelTree &geoModelTree,
    const Options &options) {
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }
  for (const auto &q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    // load data from database according to querie (Muon)
    auto qFPV = geoModelTree.geoReader
                    ->getPublishedNodes<std::string, GeoFullPhysVol *>(q);

    /** Full phyiscal volumes represent  logical detector units.*/
    for (const auto &[name, fpv] : qFPV) {
      PVConstLink physVol{fpv};
      ACTS_INFO("Convert volume "<<name);
      convertFpv(name, fpv, cache, gctx);
    }
  }
}

void GeoModelDetectorObjectFactory::convertSensitive(
    const PVConstLink &geoPV, const Transform3 &transform,
    SurfaceBoundFactory& boundFactory,
    std::vector<GeoModelSensitiveSurface> &sensitives) {
  const GeoLogVol *logVol = geoPV->getLogVol();
  const GeoShape *shape = logVol->getShape();
  int shapeId = shape->typeID();
  const std::string& name = logVol->getName();
  std::shared_ptr<const IGeoShapeConverter> converter =
      geoShapesConverters(shapeId);
  if (converter == nullptr) {
    
    throw std::runtime_error("The converter for " + printGeoShape(shape) +
                             " is a nullptr");
  }

  auto converted = converter->toSensitiveSurface(geoPV, transform, boundFactory);
  if (converted.ok()) {
    sensitives.push_back(converted.value());
    const auto &[el, sf] = converted.value();

    if (!el || !sf) {
      throw std::runtime_error(
          "The Detector Element or the Surface is a nullptr");
    }
    ACTS_VERBOSE("(successfully converted: "
                 << name << " / " << printGeoShape(shape) << " / "
                 << logVol->getMaterial()->getName() << ")");
    return;
  }
  ACTS_ERROR(name << " / " << printGeoShape(shape)
                  << " could not be converted by any converter");
}

std::vector<GeoChildNodeWithTrf>
GeoModelDetectorObjectFactory::findAllSubVolumes(const PVConstLink &vol) const {
  /// Fetch the direct children of the volume
  std::vector<GeoChildNodeWithTrf> subvolumes = getChildrenWithRef(vol, false);
  std::vector<GeoChildNodeWithTrf> sensitives;
  sensitives.reserve(subvolumes.size());
  for (auto& subvolume : subvolumes) {
    /// Check whether material or GeoNameTag satisfy the user defined patterns
    if (matches(subvolume.nodeName, subvolume.volume)) {
      sensitives.push_back(subvolume);
    }
    /// If the volume has no children nothing can be done further
    if (subvolume.volume->getNChildVols() == 0){
       continue;
    }
    /// Enter the next recursion level to check whether there're sensitive children
    std::vector<GeoChildNodeWithTrf> senssubsubvolumes =
        findAllSubVolumes(subvolume.volume);
    /* Append the found volumes to the output, but  update the transforms
     * of the nodes before. They're expressed with respect to their parent,
     * which is the grand child of the volume passed to the function call
     * -> apply on each grand child also the transform of the child. */
    std::transform(std::make_move_iterator(senssubsubvolumes.begin()),
                   std::make_move_iterator(senssubsubvolumes.end()),
                   std::back_inserter(sensitives),
                   [&subvolume](GeoChildNodeWithTrf &&volume) {
                     volume.transform = subvolume.transform * volume.transform;
                     return volume;
                   });
  }
  return sensitives;
}

bool GeoModelDetectorObjectFactory::convertBox(const std::string& name) {
  auto convB = std::ranges::any_of(m_cfg.convertBox, [&](const auto &n) {
    return name.find(n) != std::string::npos;
  });
  return convB;
}

void GeoModelDetectorObjectFactory::convertFpv(const std::string &name, 
                                               const GeoFullPhysVol *fpv, 
                                               Cache &cache,
                                               const GeometryContext &gctx) {
  const auto prevSize = cache.sensitiveSurfaces.size();
  PVConstLink physVol{fpv};

  /** Search all subvolumes that may be converted to sensitive surfaces */
  std::vector<GeoChildNodeWithTrf> surfaces = findAllSubVolumes(physVol);
  std::vector<GeoModelSensitiveSurface> sensitives;
  sensitives.reserve(surfaces.size());

  for (const auto &surface : surfaces) {
    /** Align the surface with the global position of the detector */
    const Transform3 transform = fpv->getAbsoluteTransform() * surface.transform;
    convertSensitive(surface.volume, transform, *cache.surfBoundFactory, sensitives);
  }
  cache.sensitiveSurfaces.insert(cache.sensitiveSurfaces.end(),
                                 std::make_move_iterator(sensitives.begin()), 
                                 std::make_move_iterator(sensitives.end()));
  // Extract the bounding box surrounding the surface
  if (convertBox(name)) {
    // get logVol for the shape of the volume
    const GeoShape *shape = fpv->getLogVol()->getShape();  // get shape
    const Transform3 &fpvtransform = fpv->getAbsoluteTransform();

    // convert bounding boxes with surfaces inside
    std::shared_ptr<Experimental::DetectorVolume> box =
        GeoModel::convertDetectorVolume(gctx, *shape, name, fpvtransform,
                                              sensitives);
    cache.boundingBoxes.push_back(box);
  }
  // If fpv has no subs and should not be converted to volume convert to surface
  else if (surfaces.empty() && matches(name, physVol)) {
    convertSensitive(fpv, fpv->getAbsoluteTransform(), *cache.surfBoundFactory, 
                     cache.sensitiveSurfaces);
  }

  // Set the corresponding database entry name to all sensitive surfaces
  for (auto i = prevSize; i < cache.sensitiveSurfaces.size(); ++i) {
    auto &[detEl, _] = cache.sensitiveSurfaces[i];
    detEl->setDatabaseEntryName(name);
    ACTS_VERBOSE("Set database name of the DetectorElement to "<<detEl->databaseEntryName());
  }
}
// function to determine if object fits query
bool GeoModelDetectorObjectFactory::matches(const std::string &name,
                                            const PVConstLink &physvol) const {
  if (m_cfg.nameList.empty() && m_cfg.materialList.empty()) {
    return true;
  }

  auto matchName = std::ranges::any_of(m_cfg.nameList, [&](const auto &n) {
    return name.find(n) != std::string::npos;
  });

  std::string matStr = physvol->getLogVol()->getMaterial()->getName();

  auto matchMaterial = std::ranges::any_of(
      m_cfg.materialList,
      [&](const auto &m) { return matStr.find(m) != std::string::npos; });

  bool match = matchMaterial && matchName;

  // for the fullphysvol we only check the name
  if (m_cfg.nameList.empty()) {
    return matchMaterial;
  }

  // if no material specified or we're looking at fpv judge by name only
  if (m_cfg.materialList.empty() || dynamic_pointer_cast<const GeoVFullPhysVol>(physvol)) {
    return matchName;
  }
  return match;
}
}  // namespace Acts
