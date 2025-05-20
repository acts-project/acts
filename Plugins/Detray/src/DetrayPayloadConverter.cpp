// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Detray/DetrayPayloadConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <detray/geometry/shapes/annulus2D.hpp>
#include <detray/geometry/shapes/concentric_cylinder2D.hpp>
#include <detray/geometry/shapes/cylinder2D.hpp>
#include <detray/geometry/shapes/rectangle2D.hpp>
#include <detray/geometry/shapes/ring2D.hpp>
#include <detray/geometry/shapes/trapezoid2D.hpp>
#include <detray/io/frontend/definitions.hpp>
#include <detray/io/frontend/payloads.hpp>

namespace Acts {

DetrayPayloadConverter::DetrayPayloadConverter(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {}

detray::io::transform_payload DetrayPayloadConverter::convertTransform(
    const Transform3& transform) {
  detray::io::transform_payload tfPayload;

  Eigen::Map<Vector3> tr{tfPayload.tr.data()};
  tr = transform.translation();

  Eigen::Map<SquareMatrix3> rot{tfPayload.rot.data()};
  rot = transform.linear();

  return tfPayload;
}

namespace {
using enum detray::io::shape_id;

detray::io::mask_payload convertBounds(const AnnulusBounds& annulus) {
  using enum detray::annulus2D::boundaries;
  using enum AnnulusBounds::BoundValues;

  detray::io::mask_payload payload;
  payload.shape = annulus2;
  payload.boundaries.resize(e_size);
  payload.boundaries.at(e_min_r) = annulus.get(eMinR);
  payload.boundaries.at(e_max_r) = annulus.get(eMaxR);
  payload.boundaries.at(e_min_phi_rel) = annulus.get(eMinPhiRel);
  payload.boundaries.at(e_max_phi_rel) = annulus.get(eMaxPhiRel);
  payload.boundaries.at(e_average_phi) = annulus.get(eAveragePhi);
  payload.boundaries.at(e_shift_x) = annulus.get(eOriginX);
  payload.boundaries.at(e_shift_y) = annulus.get(eOriginY);

  return payload;
}

detray::io::mask_payload convertBounds(const RectangleBounds& rectangle) {
  using enum RectangleBounds::BoundValues;
  using enum detray::rectangle2D::boundaries;

  detray::io::mask_payload payload;
  payload.shape = rectangle2;
  payload.boundaries.resize(e_size);

  double minX = rectangle.get(eMinX);
  double maxX = rectangle.get(eMaxX);
  double minY = rectangle.get(eMinY);
  double maxY = rectangle.get(eMaxY);

  if (minX != -maxX || minY != -maxY) {
    throw std::runtime_error(
        "Rectangle bounds are not symmetric, detray cannot handle this");
  }

  payload.boundaries.at(e_half_x) = maxX;
  payload.boundaries.at(e_half_y) = maxY;

  return payload;
}

detray::io::mask_payload convertBounds(const CylinderBounds& cylinder,
                                       detray::io::shape_id shape) {
  using enum CylinderBounds::BoundValues;

  detray::io::mask_payload payload;
  payload.shape = shape;
  if (shape == portal_cylinder2) {
    using enum detray::concentric_cylinder2D::boundaries;
    payload.boundaries.resize(e_size);
    payload.boundaries.at(e_r) = cylinder.get(eR);
    double hlZ = cylinder.get(eHalfLengthZ);
    payload.boundaries.at(e_lower_z) = -hlZ;
    payload.boundaries.at(e_upper_z) = hlZ;
  } else {
    using enum detray::cylinder2D::boundaries;
    payload.boundaries.resize(e_size);
    payload.boundaries.at(e_r) = cylinder.get(eR);
    double hlZ = cylinder.get(eHalfLengthZ);
    payload.boundaries.at(e_lower_z) = -hlZ;
    payload.boundaries.at(e_upper_z) = hlZ;
  }
  return payload;
}

detray::io::mask_payload convertBounds(const TrapezoidBounds& trapezoid) {
  using enum TrapezoidBounds::BoundValues;
  using enum detray::trapezoid2D::boundaries;

  detray::io::mask_payload payload;
  payload.shape = trapezoid2;

  payload.boundaries.resize(e_size);
  payload.boundaries.at(e_half_length_0) = trapezoid.get(eHalfLengthXnegY);
  payload.boundaries.at(e_half_length_1) = trapezoid.get(eHalfLengthXposY);
  payload.boundaries.at(e_half_length_2) = trapezoid.get(eHalfLengthY);
  payload.boundaries.at(e_divisor) = 1 / (2 * trapezoid.get(eHalfLengthY));

  return payload;
}

detray::io::mask_payload convertBounds(const RadialBounds& radial) {
  using enum RadialBounds::BoundValues;
  using enum detray::ring2D::boundaries;

  if (!radial.coversFullAzimuth()) {
    throw std::runtime_error(
        "Radial bounds do not cover full azimuth, detray cannot handle this");
  }

  if (radial.get(eAveragePhi) != 0.) {
    throw std::runtime_error(
        "Radial bounds have an average phi, detray cannot handle this");
  }

  detray::io::mask_payload payload;
  payload.shape = ring2;

  payload.boundaries.resize(e_size);
  payload.boundaries.at(e_inner_r) = radial.get(eMinR);
  payload.boundaries.at(e_outer_r) = radial.get(eMaxR);

  return payload;
}

}  // namespace

detray::io::mask_payload DetrayPayloadConverter::convertMask(
    const SurfaceBounds& bounds, bool forPortal) {
  detray::io::mask_payload payload;

  switch (bounds.type()) {
    using enum SurfaceBounds::BoundsType;
    using enum detray::io::shape_id;
    case eAnnulus:
      payload = convertBounds(dynamic_cast<const AnnulusBounds&>(bounds));
      break;
    case eRectangle:
      payload = convertBounds(dynamic_cast<const RectangleBounds&>(bounds));
      break;
    case eCylinder:
      payload = convertBounds(dynamic_cast<const CylinderBounds&>(bounds),
                              forPortal ? portal_cylinder2 : cylinder2);
      break;
    case eTrapezoid:
      payload = convertBounds(dynamic_cast<const TrapezoidBounds&>(bounds));
      break;
    case eDisc:
      if (auto* radial = dynamic_cast<const RadialBounds*>(&bounds);
          radial != nullptr) {
        payload = convertBounds(*radial);
      } else {
        throw std::runtime_error(
            "Disc bounds type but not radial bounds currently unsupported");
      }
      break;
    default:
      payload.shape = unknown;
      break;
  }

  return payload;
}

detray::io::surface_payload DetrayPayloadConverter::convertSurface(
    const GeometryContext& gctx, const Surface& surface, bool portal) const {
  detray::io::surface_payload payload;

  payload.transform =
      DetrayPayloadConverter::convertTransform(surface.transform(gctx));
  payload.source = surface.geometryId().value();
  payload.barcode = std::nullopt;

  bool isSensitive = false;
  if (m_cfg.sensitiveStrategy ==
      DetrayPayloadConverter::Config::SensitiveStrategy::Identifier) {
    isSensitive = surface.geometryId().sensitive() > 0;
  } else {
    isSensitive = surface.associatedDetectorElement() != nullptr;
  }

  if (portal) {
    payload.type = detray::surface_id::e_portal;
  } else {
    payload.type = isSensitive ? detray::surface_id::e_sensitive
                               : detray::surface_id::e_passive;
  }
  payload.mask = convertMask(surface.bounds(), portal);
  return payload;
}

detray::io::volume_payload DetrayPayloadConverter::convertVolume(
    const TrackingVolume& volume) const {
  detray::io::volume_payload payload;
  payload.transform =
      DetrayPayloadConverter::convertTransform(volume.transform());
  payload.name = volume.volumeName();
  switch (volume.volumeBounds().type()) {
    using enum VolumeBounds::BoundsType;
    using enum detray::volume_id;
    case eCylinder:
      payload.type = e_cylinder;
      break;
    case eCuboid:
      payload.type = e_cuboid;
      break;
    case eTrapezoid:
      payload.type = e_trapezoid;
      break;
    case eCone:
      payload.type = e_cone;
      break;
    default:
      payload.type = e_unknown;
      break;
  }
  return payload;
}

void DetrayPayloadConverter::handlePortalLink(
    const GeometryContext& gctx, const TrackingVolume& volume,
    detray::io::volume_payload& volPayload,
    std ::function<std::size_t(const TrackingVolume*)> volumeLookup,
    const PortalLinkBase& link) const {
  auto handle = [&](const TrivialPortalLink& trivial) {
    ACTS_VERBOSE("Converting trivial portal link registered to volume "
                 << volume.volumeName());
    if (&trivial.volume() == &volume) {
      ACTS_VERBOSE("~> points at this volume (" << volume.volumeName()
                                                << ") => skpping");
      return;
    }

    ACTS_VERBOSE("~> points at different volume ("
                 << trivial.volume().volumeName()
                 << ") => adding link to this volume (" << volume.volumeName()
                 << ")");

    // add the surface (including mask first)
    auto& srfPayload = volPayload.surfaces.emplace_back(
        convertSurface(gctx, trivial.surface(), true));
    srfPayload.index_in_coll = volPayload.surfaces.size() - 1;

    // lookup the target volume index (we already converted this)
    ACTS_VERBOSE("Target volume index for " << trivial.volume().volumeName()
                                            << ": "
                                            << volumeLookup(&trivial.volume()));
    auto targetVolumeIndex = volumeLookup(&trivial.volume());
    srfPayload.mask.volume_link.link = targetVolumeIndex;
  };

  if (auto* trivial = dynamic_cast<const TrivialPortalLink*>(&link);
      trivial != nullptr) {
    handle(*trivial);
  } else if (auto* composite = dynamic_cast<const CompositePortalLink*>(&link);
             composite != nullptr) {
    ACTS_VERBOSE("Converting composite portal link with "
                 << composite->links().size() << " sub-links");
    for (const auto& subLink : composite->links()) {
      const auto* subTrivial = dynamic_cast<const TrivialPortalLink*>(&subLink);

      if (subTrivial == nullptr) {
        throw std::runtime_error(
            "Composite portal link contains non-trivial portal links");
      } else {
        handle(*subTrivial);
      }
    }
  } else if (auto* grid = dynamic_cast<const GridPortalLink*>(&link);
             grid != nullptr) {
    ACTS_VERBOSE("Converting grid portal link with "
                 << grid->artifactPortalLinks().size() << " link artifacts");
    for (const auto& artifact : grid->artifactPortalLinks()) {
      handle(artifact);
    }
  } else {
    throw std::runtime_error(
        "Unknown portal link type, detray cannot handle this");
  }
}

void DetrayPayloadConverter::makeEndOfWorld(
    const GeometryContext& gctx, detray::io::volume_payload& volPayload,
    const Surface& surface) const {
  ACTS_VERBOSE("Adding end of world surface");
  auto& srfPayload =
      volPayload.surfaces.emplace_back(convertSurface(gctx, surface, true));
  srfPayload.index_in_coll = volPayload.surfaces.size() - 1;

  // Marker for end of world is MAX
  srfPayload.mask.volume_link.link = std::numeric_limits<std::size_t>::max();
}

detray::io::detector_payload DetrayPayloadConverter::convertTrackingGeometry(
    const GeometryContext& gctx, const TrackingGeometry& geometry) const {
  ACTS_INFO("Converting tracking geometry to detray format");

  if (geometry.geometryVersion() != TrackingGeometry::GeometryVersion::Gen3) {
    ACTS_WARNING(
        "Only Gen3 tracking geometries are supported. Gen1 geometries will "
        "give wrong results");
  }

  detray::io::detector_payload payload;
  std::unordered_map<const TrackingVolume*, std::size_t> volumeIds;

  auto lookup = [&volumeIds](const TrackingVolume* v) {
    return volumeIds.at(v);
  };

  geometry.apply([&](const TrackingVolume& volume) {
    auto& volPayload = payload.volumes.emplace_back(convertVolume(volume));
    volPayload.index.link = payload.volumes.size() - 1;
    volumeIds[&volume] = volPayload.index.link;

    for (auto& surface : volume.surfaces()) {
      auto& srfPayload =
          volPayload.surfaces.emplace_back(convertSurface(gctx, surface));
      srfPayload.index_in_coll = volPayload.surfaces.size() - 1;
      srfPayload.mask.volume_link.link = volPayload.index.link;
    }
  });

  // Run again over volumes, can lookup volume index from pointer now
  geometry.apply([&](const TrackingVolume& volume) {
    auto& volPayload = payload.volumes.at(volumeIds.at(&volume));

    for (const auto& portal : volume.portals()) {
      // hard-requirement: all portal links must decompose to trivial portal
      // links!

      auto* lAlong = portal.getLink(Direction::AlongNormal());
      auto* lOpposite = portal.getLink(Direction::OppositeNormal());

      if (lAlong == nullptr && lOpposite == nullptr) {
        // Sanity check: this shouldn't happen
        throw std::runtime_error("Portal link is not symmetric");
      }

      if (lAlong != nullptr) {
        handlePortalLink(gctx, volume, volPayload, lookup, *lAlong);
      } else {
        // can't both be nullptr
        assert(lOpposite != nullptr);
        makeEndOfWorld(gctx, volPayload, lOpposite->surface());
      }

      if (lOpposite != nullptr) {
        handlePortalLink(gctx, volume, volPayload, lookup, *lOpposite);
      } else {
        // can't both be nullptr
        assert(lAlong != nullptr);
        makeEndOfWorld(gctx, volPayload, lAlong->surface());
      }

      for (auto dir : {Direction::AlongNormal(), Direction::OppositeNormal()}) {
        const auto* link = portal.getLink(dir);

        if (link != nullptr) {
          handlePortalLink(gctx, volume, volPayload, lookup, *link);
        }
      }
    }

    ACTS_DEBUG("Volume " << volume.volumeName() << " has "
                         << volPayload.surfaces.size() << " surfaces");
    std::size_t nPortals =
        std::ranges::count_if(volPayload.surfaces, [](const auto& srfPayload) {
          return srfPayload.mask.volume_link.link !=
                 std::numeric_limits<std::size_t>::max();
        });
    ACTS_DEBUG("-> portals:        " << nPortals);
    ACTS_DEBUG("-> other surfaces: " << volPayload.surfaces.size() - nPortals);
  });

  ACTS_DEBUG("Collected " << payload.volumes.size() << " volumes");

  return payload;
}

}  // namespace Acts
