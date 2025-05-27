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
#include "Acts/Geometry/DetrayFwd.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
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

namespace {
std::vector<const TrivialPortalLink*> decomposeToTrivials(
    const PortalLinkBase& link, const Logger& logger) {
  std::vector<const TrivialPortalLink*> trivials;

  if (auto* trivial = dynamic_cast<const TrivialPortalLink*>(&link);
      trivial != nullptr) {
    trivials.push_back(trivial);
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
        trivials.push_back(subTrivial);
      }
    }
  } else if (auto* grid = dynamic_cast<const GridPortalLink*>(&link);
             grid != nullptr) {
    ACTS_VERBOSE("Converting grid portal link with "
                 << grid->artifactPortalLinks().size() << " link artifacts");
    for (const auto& artifact : grid->artifactPortalLinks()) {
      trivials.push_back(&artifact);
    }
  } else {
    throw std::runtime_error(
        "Unknown portal link type, detray cannot handle this");
  }

  return trivials;
}
}  // namespace

void DetrayPayloadConverter::handlePortalLink(
    const GeometryContext& gctx, const TrackingVolume& volume,
    detray::io::volume_payload& volPayload,
    const std::function<std::size_t(const TrackingVolume*)>& volumeLookup,
    std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    const PortalLinkBase& link) const {
  std::vector<const TrivialPortalLink*> trivials =
      decomposeToTrivials(link, logger());

  // If ANY of the trivials point at the current volume, we don't handle this
  // portal link at all, otherwise we would get a self-referencing volume link

  if (std::ranges::any_of(
          trivials, [&](const auto* t) { return &t->volume() == &volume; })) {
    ACTS_VERBOSE("At least one trivial link points at this volume ("
                 << volume.volumeName() << ") => skipping");
    return;
  }

  for (const auto* trivial : trivials) {
    ACTS_VERBOSE("Converting trivial portal link registered to volume "
                 << volume.volumeName());
    ACTS_VERBOSE(
        "Portal link surface is: " << trivial->surface().toStream(gctx));
    if (&trivial->volume() == &volume) {
      ACTS_VERBOSE("~> points at this volume (" << volume.volumeName()
                                                << ") => skipping");
      return;
    }

    ACTS_VERBOSE("~> points at different volume ("
                 << trivial->volume().volumeName()
                 << ") => adding link to this volume (" << volume.volumeName()
                 << ")");

    // add the surface (including mask first)
    auto& srfPayload = volPayload.surfaces.emplace_back(
        convertSurface(gctx, trivial->surface(), true));
    srfPayload.index_in_coll = volPayload.surfaces.size() - 1;

    // lookup the target volume index (we already converted this)
    ACTS_VERBOSE("Target volume index for "
                 << trivial->volume().volumeName() << ": "
                 << volumeLookup(&trivial->volume()));
    std::size_t targetVolumeIndex = volumeLookup(&trivial->volume());
    srfPayload.mask.volume_link.link = targetVolumeIndex;
    surfaceIndices[&trivial->surface()] = srfPayload.index_in_coll.value();
  }
}

void DetrayPayloadConverter::makeEndOfWorld(
    const GeometryContext& gctx, detray::io::volume_payload& volPayload,
    std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    const Surface& surface) const {
  ACTS_VERBOSE("Adding end of world surface");
  auto& srfPayload =
      volPayload.surfaces.emplace_back(convertSurface(gctx, surface, true));
  srfPayload.index_in_coll = volPayload.surfaces.size() - 1;

  // Marker for end of world is MAX
  srfPayload.mask.volume_link.link = std::numeric_limits<std::size_t>::max();

  surfaceIndices[&surface] = srfPayload.index_in_coll.value();
}

void DetrayPayloadConverter::handlePortal(
    const GeometryContext& gctx, const TrackingVolume& volume,
    detray::io::volume_payload& volPayload,
    const std::function<std::size_t(const TrackingVolume*)>& volumeLookup,
    std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    const Portal& portal) const {
  auto* lAlong = portal.getLink(Direction::AlongNormal());
  auto* lOpposite = portal.getLink(Direction::OppositeNormal());

  if (lAlong == nullptr && lOpposite == nullptr) {
    // Sanity check: this shouldn't happen
    throw std::runtime_error("Portal link is not symmetric");
  }

  if (lAlong != nullptr) {
    handlePortalLink(gctx, volume, volPayload, volumeLookup, surfaceIndices,
                     *lAlong);
  } else {
    // can't both be nullptr
    assert(lOpposite != nullptr);
    makeEndOfWorld(gctx, volPayload, surfaceIndices, lOpposite->surface());
  }

  if (lOpposite != nullptr) {
    handlePortalLink(gctx, volume, volPayload, volumeLookup, surfaceIndices,
                     *lOpposite);
  } else {
    // can't both be nullptr
    assert(lAlong != nullptr);
    makeEndOfWorld(gctx, volPayload, surfaceIndices, lAlong->surface());
  }
}

namespace {
std::optional<std::size_t> findSurfaceInVolume(
    const detray::io::volume_payload& volPayload, const Surface& surface) {
  auto srfIt =
      std::ranges::find_if(volPayload.surfaces, [&](const auto& srfPayload) {
        return srfPayload.source == surface.geometryId().value();
      });

  if (srfIt == volPayload.surfaces.end()) {
    return std::nullopt;
  }

  return std::distance(volPayload.surfaces.begin(), srfIt);
}

constexpr static detray::io::material_slab_payload s_dummyMaterialSlab{
    .type = detray::io::material_id::slab,
    .index_in_coll = std::numeric_limits<std::size_t>::max(),
    .thickness = 42,
    .mat = {42, 42, 42, 42, 42, 42, 42},
};

}  // namespace

std::pair<std::vector<detray::io::grid_payload<
              detray::io::material_slab_payload, detray::io::material_id>>,
          detray::io::material_volume_payload>
DetrayPayloadConverter::convertMaterial(
    const TrackingVolume& volume,
    const std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    detray::io::volume_payload& volPayload) const {
  ACTS_DEBUG("Converting material for volume " << volume.volumeName());
  std::vector<detray::io::grid_payload<detray::io::material_slab_payload,
                                       detray::io::material_id>>
      grids;
  detray::io::material_volume_payload homogeneous;
  homogeneous.volume_link.link = volPayload.index.link;

  ACTS_WARNING(
      "Adding dummy material slabs to homogeneous collection (detray "
      "hack)");
  for (const auto& surface : volPayload.surfaces) {
    auto& slabPayload = homogeneous.mat_slabs.emplace_back(s_dummyMaterialSlab);
    slabPayload.index_in_coll = homogeneous.mat_slabs.size() - 1;
    slabPayload.surface.link = surface.index_in_coll.value();
  }

  auto assignMaterial = [&](DetraySurfaceMaterial& detrayMaterial,
                            std::size_t srfIdx) {
    auto handleHomogeneous =
        [&](const detray::io::material_slab_payload& slab) {
          // ACTS_DEBUG("Assigning homogeneous material slab to surface "
          //            << srfIdx);
          homogeneous.mat_slabs.emplace_back(slab);
          homogeneous.mat_slabs.back().index_in_coll =
              homogeneous.mat_slabs.size() - 1;
          homogeneous.mat_slabs.back().surface.link = srfIdx;
        };

    auto handleGrid =
        [&](const detray::io::grid_payload<detray::io::material_slab_payload,
                                           detray::io::material_id>& grid) {
          ACTS_DEBUG("Assigning grid material to surface " << srfIdx);
          grids.emplace_back(grid);
          grids.back().owner_link.link = srfIdx;
        };

    std::visit(overloaded{handleHomogeneous, handleGrid}, detrayMaterial);
  };

  auto printSurfaceInfo = [&](DetraySurfaceMaterial& detrayMaterial,
                              const Surface& surface) {
    auto handleHomogeneous = [&](const detray::io::material_slab_payload&) {
      ACTS_VERBOSE("Surface " << surface.geometryId()
                              << " has homogeneous material");
    };

    auto handleGrid =
        [&](const detray::io::grid_payload<detray::io::material_slab_payload,
                                           detray::io::material_id>&) {
          ACTS_VERBOSE("Surface " << surface.geometryId()
                                  << " has grid material");
        };
    std::visit(overloaded{handleHomogeneous, handleGrid}, detrayMaterial);
  };

  for (const auto& surface : volume.surfaces()) {
    auto srfIt = surfaceIndices.find(&surface);

    if (srfIt == surfaceIndices.end()) {
      ACTS_ERROR("Surface " << surface.geometryId().value()
                            << " not found in volume " << volPayload.name
                            << ". This is a bug in the conversion.");
      throw std::runtime_error("Surface not found in volume");
    }

    std::size_t srfIdx = srfIt->second;

    if (surface.surfaceMaterial() == nullptr) {
      continue;
    }

    auto detrayMaterial = surface.surfaceMaterial()->toDetrayPayload();

    if (detrayMaterial == nullptr) {
      continue;
    }

    printSurfaceInfo(*detrayMaterial, surface);

    assignMaterial(*detrayMaterial, srfIdx);
  }

  // Portals need special treatment: we have decomposed them to their trivial
  // portal links, and only registered a subset to them as well. We again need
  // to decompose here, and only look for the portals that are actually found in
  // the volume payload
  for (const auto& portal : volume.portals()) {
    // First check, if the combined portal surface has material assigned at all,
    // if not there's nothing to do
    if (portal.surface().surfaceMaterial() == nullptr) {
      continue;
    }

    auto detrayMaterial = portal.surface().surfaceMaterial()->toDetrayPayload();

    // Portal surface material reports it does not apply to detray, skip
    if (detrayMaterial == nullptr) {
      continue;
    }

    printSurfaceInfo(*detrayMaterial, portal.surface());

    // Have valid detray material, now we need to find the surfaces that are
    // actually there in detray

    for (auto dir : {Direction::AlongNormal(), Direction::OppositeNormal()}) {
      const auto* link = portal.getLink(dir);

      if (link == nullptr) {
        continue;
      }

      std::vector<const TrivialPortalLink*> trivials =
          decomposeToTrivials(*link, logger());

      // ACTS_DEBUG("Portal link in volume " << volPayload.name << " has
      // produced "
      //                                     << trivials.size() << " trivials");
      for (const auto* trivial : trivials) {
        auto srfIt = surfaceIndices.find(&trivial->surface());

        if (srfIt == surfaceIndices.end()) {
          // This trivial was not converted, skip
          continue;
        }

        std::size_t srfIdx = srfIt->second;

        // ACTS_DEBUG("Trivial portal link in volume "
        //            << volPayload.name << " has surface "
        //            << trivial->surface().geometryId() << " detray idx "
        //            << srfIdx);

        // Assign (a copy of) the detray material to the surface payload
        // associated with this trivial
        assignMaterial(*detrayMaterial, srfIdx);
      }
    }
  }

  return {grids, homogeneous};
}

DetrayPayloadConverter::Payloads
DetrayPayloadConverter::convertTrackingGeometry(
    const GeometryContext& gctx, const TrackingGeometry& geometry) const {
  ACTS_INFO("Converting tracking geometry to detray format");

  if (geometry.geometryVersion() != TrackingGeometry::GeometryVersion::Gen3) {
    ACTS_WARNING(
        "Only Gen3 tracking geometries are supported. Gen1 geometries will "
        "give wrong results");
  }

  if (m_cfg.beampipeVolume == nullptr) {
    throw std::runtime_error("Beampipe volume not set");
  }

  Payloads payloads;
  payloads.detector = std::make_unique<detray::io::detector_payload>();
  payloads.homogeneousMaterial =
      std::make_unique<detray::io::detector_homogeneous_material_payload>();
  payloads.materialGrids = std::make_unique<detray::io::detector_grids_payload<
      detray::io::material_slab_payload, detray::io::material_id>>();

  detray::io::detector_payload& detPayload = *payloads.detector;
  detray::io::detector_homogeneous_material_payload& dthmPayload =
      *payloads.homogeneousMaterial;
  detray::io::detector_grids_payload<detray::io::material_slab_payload,
                                     detray::io::material_id>& materialGrids =
      *payloads.materialGrids;

  std::unordered_map<const TrackingVolume*, std::size_t> volumeIds;

  auto lookup = [&volumeIds](const TrackingVolume* v) {
    return volumeIds.at(v);
  };

  std::unordered_map<const TrackingVolume*,
                     std::unordered_map<const Surface*, std::size_t>>
      volumeSurfaceIndices;

  geometry.apply([&](const TrackingVolume& volume) {
    auto& volPayload = detPayload.volumes.emplace_back(convertVolume(volume));
    volPayload.index.link = detPayload.volumes.size() - 1;
    volumeIds[&volume] = volPayload.index.link;

    ACTS_DEBUG("Volume " << volume.volumeName() << " has index "
                         << volPayload.index.link);

    auto& surfaceIndices = volumeSurfaceIndices[&volume];

    for (auto& surface : volume.surfaces()) {
      auto& srfPayload =
          volPayload.surfaces.emplace_back(convertSurface(gctx, surface));
      srfPayload.index_in_coll = volPayload.surfaces.size() - 1;
      srfPayload.mask.volume_link.link = volPayload.index.link;
      surfaceIndices[&surface] = srfPayload.index_in_coll.value();
    }
  });

  // Run again over volumes, can lookup volume index from pointer now
  geometry.apply([&](const TrackingVolume& volume) {
    auto& volPayload = detPayload.volumes.at(volumeIds.at(&volume));
    auto& surfaceIndices = volumeSurfaceIndices[&volume];

    for (const auto& portal : volume.portals()) {
      handlePortal(gctx, volume, volPayload, lookup, surfaceIndices, portal);
    }

    ACTS_DEBUG("Volume " << volume.volumeName()
                         << " (idx: " << volPayload.index.link << ") has "
                         << volPayload.surfaces.size() << " total surfaces");

    std::size_t nPortals =
        std::ranges::count_if(volPayload.surfaces, [](const auto& srfPayload) {
          return srfPayload.type == detray::surface_id::e_portal;
        });
    ACTS_DEBUG("-> portals:        " << nPortals);
    std::size_t nSensitives =
        std::ranges::count_if(volPayload.surfaces, [](const auto& srfPayload) {
          return srfPayload.type == detray::surface_id::e_sensitive;
        });
    ACTS_DEBUG("-> sensitives:     " << nSensitives);
    ACTS_DEBUG("-> other surfaces: " << volPayload.surfaces.size() - nPortals -
                                            nSensitives);

    for (const auto& [surface, idx] : surfaceIndices) {
      ACTS_VERBOSE("Surface " << surface->geometryId() << " (&: " << surface
                              << ") has index " << idx);
    }

    // Portals have produced surfaces and are added in volume payload, handle
    // material now

    auto [grids, homogeneous] =
        convertMaterial(volume, surfaceIndices, volPayload);

    ACTS_DEBUG("Volume " << volume.volumeName()
                         << " (idx: " << volPayload.index.link << ") has "
                         << homogeneous.mat_slabs.size() << " material slabs");

    if (!homogeneous.mat_slabs.empty()) {
      // Only add if it's not empty (it might be)
      // NOTE: Currently, it'll always be populated by at least the homogeneous
      // NOTE: Volume association is internal to
      // `detray::io::material_volume_payload`
      dthmPayload.volumes.emplace_back(std::move(homogeneous));
    }

    ACTS_DEBUG("Volume " << volume.volumeName()
                         << " (idx: " << volPayload.index.link << ") has "
                         << grids.size() << " grids");
    if (!grids.empty()) {
      // Only add if we have grids
      // NOTE: Volume association is EXTERNAL, i.e. we need to fill a map keyed
      // by the volume index
      materialGrids.grids[volPayload.index.link] = std::move(grids);
    }
  });

  // HACK: Beampipe MUST have index 0
  // Find beampipe volume by name
  std::size_t beampipeIdx = volumeIds.at(m_cfg.beampipeVolume);
  ACTS_DEBUG("Beampipe volume (" << m_cfg.beampipeVolume->volumeName()
                                 << ") index: " << beampipeIdx);
  ACTS_DEBUG("Volume at index 0 is " << detPayload.volumes.at(0).name);

  // Swap beampipe volume to index 0
  std::swap(detPayload.volumes.at(0), detPayload.volumes.at(beampipeIdx));

  for (auto& vol : detPayload.volumes) {
    for (auto& srf : vol.surfaces) {
      if (srf.mask.volume_link.link == beampipeIdx) {
        srf.mask.volume_link.link = 0;
      } else if (srf.mask.volume_link.link == 0) {
        srf.mask.volume_link.link = beampipeIdx;
      }
    }
  }

  for (auto& mat : dthmPayload.volumes) {
    if (mat.volume_link.link == beampipeIdx) {
      mat.volume_link.link = 0;
    } else if (mat.volume_link.link == 0) {
      mat.volume_link.link = beampipeIdx;
    }
  }

  auto beampipeGridIt = materialGrids.grids.find(beampipeIdx);
  auto worldGridIt = materialGrids.grids.find(0);

  if (beampipeGridIt != materialGrids.grids.end() &&
      worldGridIt != materialGrids.grids.end()) {
    // BOTH world and beampipe have grid specifiers: swap them
    ACTS_DEBUG("Swapping beampipe and world grid specifiers");
    // auto beampipeGrid = std::move(beampipeGridIt->second);
    // materialGrids.grids.erase(beampipeGridIt);
    // auto worldGrid = std::move(worldGridIt->second);
    // materialGrids.grids.erase(worldGridIt);
    // materialGrids.grids[0] = std::move(beampipeGrid);
    // materialGrids.grids[beampipeIdx] = std::move(worldGrid);
    std::swap(beampipeGridIt->second, worldGridIt->second);
  } else if (beampipeGridIt != materialGrids.grids.end()) {
    // ONLY beampipe has grid specifier: move it to world
    ACTS_DEBUG("Moving beampipe grid specifier to world");
    materialGrids.grids[0] = std::move(beampipeGridIt->second);
    materialGrids.grids.erase(beampipeGridIt);
  } else if (worldGridIt != materialGrids.grids.end()) {
    // ONLY world has grid specifier: move it to beampipe
    ACTS_DEBUG("Moving world grid specifier to beampipe");
    materialGrids.grids[beampipeIdx] = std::move(worldGridIt->second);
    materialGrids.grids.erase(worldGridIt);
  }

  ACTS_DEBUG("Collected " << detPayload.volumes.size() << " volumes");

  return payloads;
}

}  // namespace Acts
