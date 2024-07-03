// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
// clang-format on

#include "Acts/Plugins/GeoModel/GeoModelDetectorSurfaceFactory.hpp"

#include "Acts/Plugins/GeoModel/GeoModelConverters.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"

#include <set>

namespace {
std::string recType(const GeoShapeShift &gshift);
std::string recType(const GeoShapeUnion &gunion);
std::string recType(const GeoShape &gshape);

std::string recType(const GeoShapeShift &gshift) {
  return "Shift[" + recType(*gshift.getOp()) + "]";
}
std::string recType(const GeoShapeUnion &gunion) {
  return "Union[" + recType(*gunion.getOpA()) + ", " +
         recType(*gunion.getOpB()) + "]";
}
std::string recType(const GeoShape &gshape) {
  if (auto ptr = dynamic_cast<const GeoShapeUnion *>(&gshape); ptr != nullptr) {
    return recType(*ptr);
  }
  if (auto ptr = dynamic_cast<const GeoShapeShift *>(&gshape); ptr != nullptr) {
    return recType(*ptr);
  }
  return gshape.type();
}

}  // namespace

Acts::GeoModelDetectorSurfaceFactory::GeoModelDetectorSurfaceFactory(
    const Config &cfg, std::unique_ptr<const Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {}

void Acts::GeoModelDetectorSurfaceFactory::construct(
    Cache &cache, const GeometryContext &, const GeoModelTree &geoModelTree,
    const Options &options) {
  if (geoModelTree.geoReader == nullptr) {
    throw std::invalid_argument("GeoModelTree has no GeoModelReader");
  }

  for (const auto &q : options.queries) {
    ACTS_VERBOSE("Constructing detector elements for query " << q);
    auto qFPV = geoModelTree.geoReader
                    ->getPublishedNodes<std::string, GeoFullPhysVol *>(q);

    auto matches = [&](const std::string &name, PVConstLink physvol) {
      if (m_cfg.nameList.empty() && m_cfg.materialList.empty()) {
        return true;
      }

      bool matchName{false}, matchMaterial{false};

      matchName |= std::any_of(
          m_cfg.nameList.begin(), m_cfg.nameList.end(),
          [&](const auto &n) { return name.find(n) != std::string::npos; });

      std::string matStr = physvol->getLogVol()->getMaterial()->getName();

      matchMaterial |= std::any_of(
          m_cfg.materialList.begin(), m_cfg.materialList.end(),
          [&](const auto &m) { return matStr.find(m) != std::string::npos; });

      bool match = matchMaterial && matchName;
      GeoIntrusivePtr<const GeoVFullPhysVol> fullVol =
          dynamic_pointer_cast<const GeoVFullPhysVol>(physvol);

      // for the fullphysvol we only check the name

      if (m_cfg.nameList.empty()) {
        return matchMaterial;
      }

      if (m_cfg.materialList.empty() || !(fullVol == nullptr)) {
        return matchName;
      }

      return match;
    };

    // Store stems of volume names and materials for INFO output
    std::set<std::string> volNameStems;
    std::set<std::string> materials;

    for (auto &[name, fpv] : qFPV) {
      PVConstLink physVol{fpv};
      GeoIntrusivePtr<const GeoVPhysVol> physvol{fpv};

      if (!matches(name, physVol)) {
        continue;
      }

	  // TODO: This should be changed in a way that both top level volume and subvolumes can be converted.
      if (m_cfg.convertSubVolumes) {
        std::vector<GeoChildNodeWithTrf> subvolumes =
            findAllSubVolumes(physVol, matches);
        ACTS_VERBOSE("Found" << subvolumes.size()
                             << "sub volumes matching the query");

        for (auto &subvol : subvolumes) {
          // get the transform of the subvolume to the global frame
          const Transform3 &transform =
              fpv->getAbsoluteTransform() * subvol.transform;

          convertSensitive(subvol.volume, transform, cache.sensitiveSurfaces);
        }

      } else {
        // convert the full phys vol
        convertSensitive(physVol, fpv->getAbsoluteTransform(nullptr),
                         cache.sensitiveSurfaces);
      }
    }

    ACTS_VERBOSE("Found " << qFPV.size()
                          << " full physical volumes matching the query.");
  }

  ACTS_DEBUG("Constructed "
             << cache.sensitiveSurfaces.size() << " sensitive elements and "
             << cache.passiveSurfaces.size() << " passive elements");
}

void Acts::GeoModelDetectorSurfaceFactory::convertSensitive(
    PVConstLink geoPV, const Acts::Transform3 &transform,
    std::vector<Acts::GeoModelSensitiveSurface> &sensitives) {
  const GeoLogVol *logVol = geoPV->getLogVol();
  const GeoShape *shape = logVol->getShape();
  int shapeId = shape->typeID();
  std::string name = logVol->getName();
  std::shared_ptr<const Acts::IGeoShapeConverter> converter =
      Acts::GeoShapesConverters(shapeId);
  if (converter != nullptr) {
    auto converted = converter->toSensitiveSurface(geoPV, transform);
    if (converted.ok()) {
      sensitives.push_back(converted.value());
      const auto &[el, sf] = converted.value();

      ACTS_VERBOSE("(successfully converted: "
                   << name << " / " << recType(*shape) << " / "
                   << logVol->getMaterial()->getName() << ")");

      if (!(el && sf)) {
        throw std::runtime_error("something is nullptr");
      }
      return;
    }
    ACTS_DEBUG(name << " / " << recType(*shape)
                    << ") could not be converted by any converter");
  }
  ACTS_DEBUG("converter is a nullptr");
  return;
}

std::vector<GeoChildNodeWithTrf>
Acts::GeoModelDetectorSurfaceFactory::findAllSubVolumes(
    PVConstLink geoPV,
    std::function<bool(std::string, PVConstLink)> matchFunc) {
  const std::vector<GeoChildNodeWithTrf> children =
      getChildrenWithRef(geoPV, false);
  std::vector<GeoChildNodeWithTrf> foundVols{};

  for (const auto &child : children) {
    if (matchFunc(child.nodeName, child.volume)) {
      ACTS_VERBOSE("The subvol" << child.nodeName << "matches the queries");
      foundVols.push_back(child);
    }
    if (child.volume->getNChildVols() != 0) {
      continue;
    }
    std::vector<GeoChildNodeWithTrf> grandChildren =
        findAllSubVolumes(child.volume, matchFunc);
    std::transform(std::make_move_iterator(grandChildren.begin()),
                   std::make_move_iterator(grandChildren.end()),
                   std::back_inserter(foundVols),
                   [&child](GeoChildNodeWithTrf &&vol) {
                     vol.transform = child.transform * vol.transform;
                     return vol;
                   });
  }

  return foundVols;
}
