// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/GeoModel/GeoModelDetectorSurfaceFactory.hpp"

#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Plugins/GeoModel/converters/GeoUnionDoubleTrdConverter.hpp"

#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoShapeUnion.h>

#include <set>

namespace {
std::string gname(const GeoShapeShift &gshift);
std::string gname(const GeoShapeUnion &gunion);
std::string gname(const GeoShape &gshape);

std::string gname(const GeoShapeShift &gshift) {
  return "Shift[" + gname(*gshift.getOp()) + "]";
}
std::string gname(const GeoShapeUnion &gunion) {
  return "Union[" + gname(*gunion.getOpA()) + ", " + gname(*gunion.getOpB()) +
         "]";
}
std::string gname(const GeoShape &gshape) {
  if (auto ptr = dynamic_cast<const GeoShapeUnion *>(&gshape); ptr != nullptr) {
    return gname(*ptr);
  }
  if (auto ptr = dynamic_cast<const GeoShapeShift *>(&gshape); ptr != nullptr) {
    return gname(*ptr);
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

    auto matches = [&](const std::string &name, const GeoVFullPhysVol &fpv) {
      bool match = m_cfg.nameList.empty() && m_cfg.materialList.empty();
      if(match){
        return true;
      }

      match |= std::any_of(
          m_cfg.nameList.begin(), m_cfg.nameList.end(),
          [&](const auto &n) { return name.find(n) != std::string::npos; });
      if(match){
        return true;
      }

      const auto &matStr = fpv.getLogVol()->getMaterial()->getName();
      match |= std::any_of(
          m_cfg.materialList.begin(), m_cfg.materialList.end(),
          [&](const auto &m) { return matStr.find(m) != std::string::npos; });

      return match;
    };

    std::set<std::string> nameStems;
    std::set<std::string> materials;

    for (auto &[name, fpv] : qFPV) {
      const std::string &vname = fpv->getLogVol()->getName();
      const GeoShape &shape = *fpv->getLogVol()->getShape();
      // Mask
      if (!matches(name, *fpv)) {
        continue;
      }

      bool success = false;

      // This is only hacked right now, for a proper solution we should
      // restructure the code and add a convert function that recursively calls
      // itself for unions
      if (auto gunion = dynamic_cast<const GeoShapeUnion *>(&shape);
          gunion != nullptr) {
        GeoUnionConverter converter;
        converter.useA = true;
        auto converted1 = converter.toSensitiveSurface(*fpv);
        if (converted1.ok()) {
          cache.sensitiveSurfaces.push_back(converted1.value());
        }
        converter.useA = false;
        auto converted2 = converter.toSensitiveSurface(*fpv);
        if (converted2.ok()) {
          cache.sensitiveSurfaces.push_back(converted2.value());
        }
        success = converted1.ok() && converted2.ok();
        if (success) {
          ACTS_VERBOSE("successfully converted " << name << " (" << vname
                                                 << " / " << gname(shape)
                                                 << ")      mat: " << fpv->getLogVol()->getMaterial()->getName());

        auto isBigTrap = detail::unionIsBigTrapezoid(*std::get<1>(converted1.value()), *std::get<1>(converted2.value()));
        ACTS_DEBUG("Is union big trapezoid? " << (isBigTrap != nullptr));
        if (isBigTrap) {
          cache.sensitiveSurfaces.push_back(GeoModelSensitiveSurface{nullptr, isBigTrap});
        }
        }
      }
      // Convert using the list of converters
      else {
        for (const auto &converter : m_cfg.shapeConverters) {
          auto converted = converter->toSensitiveSurface(*fpv);
          if (converted.ok()) {
            // Add the element and surface to the cache
            cache.sensitiveSurfaces.push_back(converted.value());
            success = true;
            ACTS_VERBOSE("successfully converted " << name << " (" << vname
                                                   << " / " << gname(shape)
                                                   << ")      mat: " << fpv->getLogVol()->getMaterial()->getName());
            break;
          }
        }
      }

      if (!success) {
        ACTS_DEBUG(name << " (" << vname << " / " << gname(shape)
                        << ") could not be converted by any converter");
      } else {
        nameStems.emplace(name.substr(0, 6));
        materials.emplace(fpv->getLogVol()->getMaterial()->getName());
      }
    }
    ACTS_VERBOSE("Found " << qFPV.size()
                          << " full physical volumes matching the query.");
    {
      ACTS_DEBUG("Stems of converted volumes: " << [&]() {
        std::stringstream ss;
        for (const auto &el : nameStems) {
          ss << el << " ";
        }
        return ss.str();
      }());
      ACTS_DEBUG("Converted volumes with materials: " << [&]() {
        std::stringstream ss;
        for (const auto &el : materials) {
          ss << el << " ";
        }
        return ss.str();
      }());
    }
  }
  ACTS_DEBUG("Constructed "
             << cache.sensitiveSurfaces.size() << " sensitive elements and "
             << cache.passiveSurfaces.size() << " passive elements");
}
