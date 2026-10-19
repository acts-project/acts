// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"

#include <algorithm>
#include <limits>
#include <optional>
#include <set>

#include <detray/geometry/shapes/annulus2D.hpp>
#include <detray/geometry/shapes/concentric_cylinder2D.hpp>
#include <detray/geometry/shapes/cylinder2D.hpp>
#include <detray/geometry/shapes/rectangle2D.hpp>
#include <detray/geometry/shapes/ring2D.hpp>
#include <detray/geometry/shapes/trapezoid2D.hpp>
#include <detray/io/frontend/definitions.hpp>
#include <detray/io/frontend/payloads.hpp>

using namespace Acts;

namespace ActsPlugins {

DetrayPayloadConverter::DetrayPayloadConverter(
    const Config& config, std::unique_ptr<const Logger> logger)
    : m_cfg(config), m_logger(std::move(logger)) {}

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

  payload.transform = DetrayConversionUtils::convertTransform(
      surface.localToGlobalTransform(gctx));
  payload.source = surface.geometryId().value();
  payload.identifier = std::nullopt;

  bool isSensitive = false;
  if (m_cfg.sensitiveStrategy ==
      DetrayPayloadConverter::Config::SensitiveStrategy::Identifier) {
    isSensitive = surface.geometryId().sensitive() > 0;
  } else {
    isSensitive = surface.isSensitive();
  }

  if (portal) {
    payload.type = detray::surface_id::e_portal;
  } else {
    payload.type = isSensitive ? detray::surface_id::e_sensitive
                               : detray::surface_id::e_passive;
  }
  payload.masks = {convertMask(surface.bounds(), portal)};
  return payload;
}

detray::io::volume_payload DetrayPayloadConverter::convertVolume(
    const GeometryContext& gctx, const TrackingVolume& volume) const {
  detray::io::volume_payload payload;
  payload.transform = DetrayConversionUtils::convertTransform(
      volume.localToGlobalTransform(gctx));
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

/// A contiguous piece of a portal along its segmentation direction, expressed
/// in the local frame of the portal surface, together with the volume the
/// portal leads into on that piece (nullptr for the end of the world).
struct PortalSegment {
  double min;
  double max;
  const TrackingVolume* volume;
};

/// Direction along which a portal is split into segments. Cylinders are split
/// in z and discs in r, which is what detray's multi-mask portals support. For
/// planes, the direction is taken from the first composite or 1D grid link.
AxisDirection segmentationDirection(const Portal& portal) {
  switch (portal.surface().type()) {
    using enum Surface::SurfaceType;
    case Cylinder:
      return AxisDirection::AxisZ;
    case Disc:
      return AxisDirection::AxisR;
    case Plane:
      for (auto dir : {Direction::AlongNormal(), Direction::OppositeNormal()}) {
        const auto* link = portal.getLink(dir);
        if (const auto* composite =
                dynamic_cast<const CompositePortalLink*>(link);
            composite != nullptr) {
          return composite->direction();
        }
        if (const auto* grid = dynamic_cast<const GridPortalLink*>(link);
            grid != nullptr && grid->dim() == 1) {
          return grid->direction();
        }
      }
      return AxisDirection::AxisX;
    default:
      throw std::runtime_error(
          "Portal surface type is not supported by the detray conversion");
  }
}

/// Offset of @p surface along @p direction, in the local frame of @p portalSurface
double offsetAlong(const GeometryContext& gctx, const Surface& portalSurface,
                   const Surface& surface, AxisDirection direction) {
  const Vector3 offset = (portalSurface.localToGlobalTransform(gctx).inverse() *
                          surface.localToGlobalTransform(gctx))
                             .translation();
  switch (direction) {
    using enum AxisDirection;
    case AxisZ:
      return offset[eZ];
    case AxisX:
      return offset[eX];
    case AxisY:
      return offset[eY];
    default:
      // Radial extents of coplanar discs do not depend on the frame
      return 0.;
  }
}

/// Extent of @p surface along @p direction, in the local frame of @p portalSurface
std::pair<double, double> extentAlong(const GeometryContext& gctx,
                                      const Surface& portalSurface,
                                      const Surface& surface,
                                      AxisDirection direction) {
  const double offset = offsetAlong(gctx, portalSurface, surface, direction);

  if (const auto* cylinder =
          dynamic_cast<const CylinderBounds*>(&surface.bounds());
      cylinder != nullptr && direction == AxisDirection::AxisZ) {
    const double hlZ = cylinder->get(CylinderBounds::eHalfLengthZ);
    return {offset - hlZ, offset + hlZ};
  }
  if (const auto* radial = dynamic_cast<const RadialBounds*>(&surface.bounds());
      radial != nullptr && direction == AxisDirection::AxisR) {
    return {radial->get(RadialBounds::eMinR), radial->get(RadialBounds::eMaxR)};
  }
  if (const auto* rectangle =
          dynamic_cast<const RectangleBounds*>(&surface.bounds());
      rectangle != nullptr) {
    if (direction == AxisDirection::AxisX) {
      return {offset + rectangle->get(RectangleBounds::eMinX),
              offset + rectangle->get(RectangleBounds::eMaxX)};
    }
    if (direction == AxisDirection::AxisY) {
      return {offset + rectangle->get(RectangleBounds::eMinY),
              offset + rectangle->get(RectangleBounds::eMaxY)};
    }
  }

  throw std::runtime_error("Cannot segment portal surface along " +
                           axisDirectionName(direction) +
                           ": unsupported surface bounds");
}

/// Collect the segments of @p link along @p direction in the frame of
/// @p portalSurface. Grid links are read from their bins, so this does not
/// depend on the trivial links a grid was originally built from.
void collectSegments(const GeometryContext& gctx, const Surface& portalSurface,
                     const PortalLinkBase& link, AxisDirection direction,
                     std::vector<PortalSegment>& segments) {
  if (const auto* trivial = dynamic_cast<const TrivialPortalLink*>(&link);
      trivial != nullptr) {
    auto [min, max] =
        extentAlong(gctx, portalSurface, trivial->surface(), direction);
    segments.push_back({min, max, &trivial->volume()});
  } else if (const auto* composite =
                 dynamic_cast<const CompositePortalLink*>(&link);
             composite != nullptr) {
    for (const auto& child : composite->links()) {
      collectSegments(gctx, portalSurface, child, direction, segments);
    }
  } else if (const auto* grid = dynamic_cast<const GridPortalLink*>(&link);
             grid != nullptr) {
    AnyGridConstView<const TrackingVolume*> view(grid->grid());
    const auto nBins = view.multiAxisAny().getNBinsAny();

    if (grid->dim() == 1 && grid->direction() == direction) {
      const double offset =
          offsetAlong(gctx, portalSurface, grid->surface(), direction);
      const std::vector<double> edges =
          grid->grid().axes().front()->getBinEdges();
      for (std::size_t i = 0; i < nBins.at(0); ++i) {
        const TrackingVolume* target = view.atLocalBins({i + 1});
        if (target != nullptr) {
          segments.push_back(
              {edges[i] + offset, edges[i + 1] + offset, target});
        }
      }
      return;
    }

    // Binning that detray portal masks cannot express is only convertible if
    // the grid leads into a single volume everywhere.
    std::set<const TrackingVolume*> targets;
    if (grid->dim() == 1) {
      for (std::size_t i0 = 1; i0 <= nBins.at(0); ++i0) {
        targets.insert(view.atLocalBins({i0}));
      }
    } else {
      for (std::size_t i0 = 1; i0 <= nBins.at(0); ++i0) {
        for (std::size_t i1 = 1; i1 <= nBins.at(1); ++i1) {
          targets.insert(view.atLocalBins({i0, i1}));
        }
      }
    }
    targets.erase(nullptr);
    if (targets.size() != 1) {
      throw std::runtime_error(
          "Grid portal link binned along " +
          axisDirectionName(grid->direction()) + " leads into " +
          std::to_string(targets.size()) +
          " volumes, detray portals can only be segmented along " +
          axisDirectionName(direction));
    }
    auto [min, max] =
        extentAlong(gctx, portalSurface, grid->surface(), direction);
    segments.push_back({min, max, *targets.begin()});
  } else {
    throw std::runtime_error(
        "Unknown portal link type, detray cannot handle this");
  }
}

/// Sort @p segments, merge neighbours leading into the same volume, and absorb
/// segments and gaps below @p tolerance into their neighbours.
void normalizeSegments(std::vector<PortalSegment>& segments, double tolerance) {
  std::ranges::sort(segments, {}, &PortalSegment::min);

  std::vector<PortalSegment> result;
  result.reserve(segments.size());
  for (const auto& segment : segments) {
    if (result.empty()) {
      result.push_back(segment);
      continue;
    }
    auto& previous = result.back();
    const bool contiguous = segment.min - previous.max <= tolerance;
    if (contiguous && segment.volume == previous.volume) {
      previous.max = std::max(previous.max, segment.max);
    } else if (contiguous && segment.max - segment.min < tolerance) {
      previous.max = std::max(previous.max, segment.max);
    } else if (contiguous && previous.max - previous.min < tolerance) {
      previous = {previous.min, segment.max, segment.volume};
    } else {
      if (contiguous) {
        // Close numerical gaps (or overlaps) between neighbours
        previous.max = segment.min;
      }
      result.push_back(segment);
    }
  }
  segments = std::move(result);
}

/// Check compatibility between ACTS surface type and detray material_id
/// @returns pair<is_compatible, expected_material_id>
std::pair<bool, detray::io::material_id> isGridMaterialCompatible(
    Surface::SurfaceType surfaceType, detray::io::material_id materialId) {
  using enum Surface::SurfaceType;
  using enum detray::io::material_id;

  switch (surfaceType) {
    case Cylinder:
      return {materialId == concentric_cylinder2_map, concentric_cylinder2_map};
    case Disc:
      return {materialId == ring2_map, ring2_map};
    case Plane:
      return {materialId == rectangle2_map, rectangle2_map};
    default:
      // For other surface types, we don't have specific checks yet
      return {true, unknown};
  }
}

}  // namespace

void DetrayPayloadConverter::handlePortal(
    const GeometryContext& gctx, const TrackingVolume& volume,
    detray::io::volume_payload& volPayload,
    const std::function<std::size_t(const TrackingVolume*)>& volumeLookup,
    std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    const Portal& portal) const {
  const auto* lAlong = portal.getLink(Direction::AlongNormal());
  const auto* lOpposite = portal.getLink(Direction::OppositeNormal());

  if (lAlong == nullptr && lOpposite == nullptr) {
    // Sanity check: this shouldn't happen
    throw std::runtime_error("Portal link is not symmetric");
  }

  const RegularSurface& portalSurface = portal.surface();
  const AxisDirection direction = segmentationDirection(portal);

  auto segmentsOf = [&](const PortalLinkBase* link) {
    std::vector<PortalSegment> segments;
    if (link != nullptr) {
      collectSegments(gctx, portalSurface, *link, direction, segments);
    }
    return segments;
  };

  std::vector<PortalSegment> along = segmentsOf(lAlong);
  std::vector<PortalSegment> opposite = segmentsOf(lOpposite);

  auto leadsHere = [&](const PortalSegment& s) { return s.volume == &volume; };

  // The link that leads into this volume tells us which part of the portal
  // borders it. The other link provides the neighbours on that part.
  const bool alongIsOwn = std::ranges::any_of(along, leadsHere);
  if (!alongIsOwn && !std::ranges::any_of(opposite, leadsHere)) {
    ACTS_ERROR("Portal on " << portalSurface.geometryId() << " of volume "
                            << volume.volumeName()
                            << " does not lead into that volume");
    throw std::runtime_error("Portal does not lead into its volume");
  }
  const std::vector<PortalSegment>& own = alongIsOwn ? along : opposite;
  const std::vector<PortalSegment>& other = alongIsOwn ? opposite : along;
  const bool endOfWorld = (alongIsOwn ? lOpposite : lAlong) == nullptr;

  // Clip every neighbour segment to the parts of the portal that border this
  // volume. This keeps detray portal masks inside the volume they belong to.
  std::vector<PortalSegment> segments;
  for (const auto& ownSegment : own) {
    if (!leadsHere(ownSegment)) {
      continue;
    }
    if (endOfWorld) {
      segments.push_back({ownSegment.min, ownSegment.max, nullptr});
      continue;
    }
    for (const auto& otherSegment : other) {
      if (leadsHere(otherSegment)) {
        // Would be a self-referencing volume link
        continue;
      }
      const double min = std::max(ownSegment.min, otherSegment.min);
      const double max = std::min(ownSegment.max, otherSegment.max);
      if (max > min) {
        segments.push_back({min, max, otherSegment.volume});
      }
    }
  }

  normalizeSegments(segments, m_cfg.portalSegmentTolerance);

  // Parts of this volume's face without a neighbour become holes in detray
  double ownLength = 0.;
  for (const auto& ownSegment : own) {
    if (leadsHere(ownSegment)) {
      ownLength += ownSegment.max - ownSegment.min;
    }
  }
  double coveredLength = 0.;
  for (const auto& segment : segments) {
    coveredLength += segment.max - segment.min;
  }
  if (ownLength - coveredLength > m_cfg.portalSegmentTolerance) {
    ACTS_WARNING("Portal on " << portalSurface.geometryId() << " of volume "
                              << volume.volumeName() << " only covers "
                              << coveredLength << " of " << ownLength
                              << " along " << axisDirectionName(direction)
                              << " with neighbours");
  }

  if (segments.empty()) {
    ACTS_VERBOSE("Portal on " << portalSurface.geometryId()
                              << " has no neighbour for volume "
                              << volume.volumeName() << " => skipping");
    return;
  }

  ACTS_VERBOSE("Portal on " << portalSurface.geometryId() << " of volume "
                            << volume.volumeName() << " is split into "
                            << segments.size() << " segment(s) along "
                            << axisDirectionName(direction));
  for (const auto& segment : segments) {
    ACTS_VERBOSE("~> [" << segment.min << ", " << segment.max << "] -> "
                        << (segment.volume != nullptr
                                ? segment.volume->volumeName()
                                : std::string{"end of world"}));
  }

  auto linkOf = [&](const PortalSegment& segment) -> std::size_t {
    return segment.volume != nullptr ? volumeLookup(segment.volume)
                                     : std::numeric_limits<std::size_t>::max();
  };

  if (direction == AxisDirection::AxisZ || direction == AxisDirection::AxisR) {
    // Concentric cylinders and rings: all segments share the portal surface
    // transform, so a single detray portal with one mask per segment is used.
    // Surface material defined on the portal surface lines up with it as is.
    auto& srfPayload = volPayload.surfaces.emplace_back(
        convertSurface(gctx, portalSurface, true));
    srfPayload.index_in_coll = volPayload.surfaces.size() - 1;

    const detray::io::mask_payload baseMask = srfPayload.masks.at(0);
    srfPayload.masks.clear();
    for (const auto& segment : segments) {
      auto& mask = srfPayload.masks.emplace_back(baseMask);
      if (direction == AxisDirection::AxisZ) {
        using enum detray::concentric_cylinder2D::boundaries;
        mask.boundaries.at(e_lower_z) = segment.min;
        mask.boundaries.at(e_upper_z) = segment.max;
      } else {
        using enum detray::ring2D::boundaries;
        mask.boundaries.at(e_inner_r) = segment.min;
        mask.boundaries.at(e_outer_r) = segment.max;
      }
      mask.volume_link.link = linkOf(segment);
    }

    surfaceIndices[&portalSurface] = srfPayload.index_in_coll.value();
    return;
  }

  // Planes: detray rectangles are centered on their surface, so every segment
  // gets its own surface, shifted along the segmentation direction.
  if (segments.size() > 1 && portalSurface.hasMaterial()) {
    ACTS_ERROR("Plane portal on " << portalSurface.geometryId()
                                  << " carries material but is split into "
                                  << segments.size() << " segments in volume "
                                  << volume.volumeName());
    throw DetrayUnsupportedMaterialException(
        "Material on segmented plane portals is not supported");
  }

  const Transform3& portalTransform =
      portalSurface.localToGlobalTransform(gctx);
  for (const auto& segment : segments) {
    const double center = 0.5 * (segment.min + segment.max);
    const double halfLength = 0.5 * (segment.max - segment.min);
    const bool dirX = direction == AxisDirection::AxisX;

    auto& srfPayload = volPayload.surfaces.emplace_back(
        convertSurface(gctx, portalSurface, true));
    srfPayload.index_in_coll = volPayload.surfaces.size() - 1;
    srfPayload.transform = DetrayConversionUtils::convertTransform(
        portalTransform *
        Translation3{dirX ? Vector3{center, 0, 0} : Vector3{0, center, 0}});

    auto& mask = srfPayload.masks.at(0);
    using enum detray::rectangle2D::boundaries;
    mask.boundaries.at(dirX ? e_half_x : e_half_y) = halfLength;
    mask.volume_link.link = linkOf(segment);

    surfaceIndices[&portalSurface] = srfPayload.index_in_coll.value();
  }
}

std::pair<std::vector<detray::io::grid_payload<
              detray::io::surface_material_payload, detray::io::material_id>>,
          detray::io::material_volume_payload>
DetrayPayloadConverter::convertMaterial(
    const TrackingVolume& volume,
    const std::unordered_map<const Surface*, std::size_t>& surfaceIndices,
    detray::io::volume_payload& volPayload) const {
  ACTS_DEBUG("Converting material for volume " << volume.volumeName());
  std::vector<detray::io::grid_payload<detray::io::surface_material_payload,
                                       detray::io::material_id>>
      grids;
  detray::io::material_volume_payload homogeneous;
  homogeneous.volume_link.link = volPayload.index.link;

  std::map<std::size_t, const ISurfaceMaterial*> srfIdxToMaterial;

  auto assignMaterial = [&](const ISurfaceMaterial* material,
                            DetraySurfaceMaterial& detrayMaterial,
                            std::size_t srfIdx, const Surface& surface) {
    auto handleHomogeneous =
        [&](const detray::io::surface_material_payload& slab) {
          // A surface can be visited more than once (e.g. through decomposed
          // portal links), so update an existing entry instead of duplicating.
          //
          // `index_in_coll` is deliberately left unset: it is the position in
          // the detector's material collection, and forcing it to the surface
          // index would make detray size that collection to the largest
          // surface index and default-fill the gaps with invalid material.
          // Leaving it unset packs the collection in payload order instead.
          auto it = std::ranges::find_if(
              homogeneous.surface_mat, [srfIdx](const auto& matslab) {
                return matslab.surface.link == srfIdx;
              });

          if (it != homogeneous.surface_mat.end()) {
            ACTS_VERBOSE("Updating slab in homogeneous material for surface "
                         << srfIdx);
            auto& targetSlab = *it;
            targetSlab = slab;
            targetSlab.index_in_coll.reset();
            targetSlab.surface.link = srfIdx;
          } else {
            ACTS_VERBOSE("Adding slab to homogeneous material for surface "
                         << srfIdx);
            auto& newSlab = homogeneous.surface_mat.emplace_back(slab);
            newSlab.index_in_coll.reset();
            newSlab.surface.link = srfIdx;
          }

          auto sit = srfIdxToMaterial.find(srfIdx);
          if (sit != srfIdxToMaterial.end() && sit->second != material) {
            ACTS_ERROR("Surface "
                       << srfIdx
                       << " already has a different material assigned");
            throw std::runtime_error("Material mismatch for surface");
          }
        };

    auto handleGrid = [&](const detray::io::grid_payload<
                          detray::io::surface_material_payload,
                          detray::io::material_id>& grid) {
      ACTS_DEBUG("Assigning grid material to surface " << srfIdx);
      auto it = srfIdxToMaterial.find(srfIdx);

      if (it != srfIdxToMaterial.end()) {
        // Surface already has some material assigned
        if (it->second != material) {
          ACTS_ERROR("Surface "
                     << srfIdx << " already has a different material assigned");
          throw std::runtime_error("Material mismatch for surface");
        }

        // It's the same material again, we just skip adding a duplicate
        return;
      }

      // @TODO: Add a consistency check between the surface type and the
      //        axis types: detray's consistency check will find this but
      //        the error message is not trivial. In this location, we can
      //        print out the surface id and guide debugging!

      // Check compatibility between surface type and grid material_id
      auto [isCompatible, expectedId] =
          isGridMaterialCompatible(surface.type(), grid.grid_link.type);

      if (!isCompatible) {
        ACTS_ERROR("Grid material compatibility error for surface "
                   << surface.geometryId() << " (index " << srfIdx << "): "
                   << "Surface type " << surface.type()
                   << " is incompatible with material_id '"
                   << toUnderlying(grid.grid_link.type) << "'. Expected '"
                   << toUnderlying(expectedId)
                   << "'. This indicates the BinUtility axis definition "
                   << "does not match the surface geometry.");
        throw std::runtime_error(
            "Grid material has incompatible axis definition for surface type");
      }

      grids.emplace_back(grid);
      grids.back().owner_link.link = srfIdx;
    };

    std::visit(overloaded{handleHomogeneous, handleGrid}, detrayMaterial);

    srfIdxToMaterial[srfIdx] = material;
  };

  auto printSurfaceInfo = [&](DetraySurfaceMaterial& detrayMaterial,
                              const Surface& surface) {
    auto handleHomogeneous = [&](const detray::io::surface_material_payload&) {
      ACTS_VERBOSE("Surface " << surface.geometryId()
                              << " has homogeneous material");
    };

    auto handleGrid =
        [&](const detray::io::grid_payload<detray::io::surface_material_payload,
                                           detray::io::material_id>&) {
          ACTS_VERBOSE("Surface " << surface.geometryId()
                                  << " has grid material");
        };
    std::visit(overloaded{handleHomogeneous, handleGrid}, detrayMaterial);
  };

  ACTS_VERBOSE("Looping over " << volume.surfaces().size()
                               << " surfaces in volume " << volPayload.name);
  for (const auto& surface : volume.surfaces()) {
    auto srfIt = surfaceIndices.find(&surface);

    if (srfIt == surfaceIndices.end()) {
      ACTS_ERROR("Surface " << surface.geometryId().value()
                            << " not found in volume " << volPayload.name
                            << ". This is a bug in the conversion.");
      throw std::runtime_error("Surface not found in volume");
    }

    std::size_t srfIdx = srfIt->second;

    if (!surface.hasMaterial()) {
      continue;
    }

    std::optional detrayMaterial =
        m_cfg.convertSurfaceMaterial(*surface.surfaceMaterial(), surface);

    if (!detrayMaterial.has_value()) {
      continue;
    }

    printSurfaceInfo(*detrayMaterial, surface);

    assignMaterial(surface.surfaceMaterial(), *detrayMaterial, srfIdx, surface);
  }
  ACTS_VERBOSE("Looping over " << volume.portals().size()
                               << " portals in volume " << volPayload.name);

  // Cylinder and disc portals are converted into a single detray surface per
  // volume, which carries the portal surface transform and one mask per
  // neighbour. The material defined on the portal surface applies to it as is.
  // Segmented plane portals with material are rejected in handlePortal.
  for (const auto& portal : volume.portals()) {
    // First check, if the portal surface has material assigned at all, if not
    // there's nothing to do
    const auto* surfaceMaterial = portal.surface().surfaceMaterial();
    if (surfaceMaterial == nullptr) {
      continue;
    }

    auto srfIt = surfaceIndices.find(&portal.surface());
    if (srfIt == surfaceIndices.end()) {
      // The portal did not produce a surface in this volume
      continue;
    }

    std::optional detrayMaterial =
        m_cfg.convertSurfaceMaterial(*surfaceMaterial, portal.surface());

    // Portal surface material reports it does not apply to detray, skip
    if (!detrayMaterial.has_value()) {
      continue;
    }

    printSurfaceInfo(*detrayMaterial, portal.surface());

    ACTS_VERBOSE("Portal on surface " << portal.surface().geometryId()
                                      << " in volume " << volPayload.name
                                      << " has detray idx " << srfIt->second);

    // Assign (a copy of) the detray material to the portal surface payload
    assignMaterial(surfaceMaterial, *detrayMaterial, srfIt->second,
                   portal.surface());
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
    throw std::runtime_error(
        "Beampipe volume not set. This is needed to ensure detray receives the "
        "beampip volume where it expects it");
  }

  Payloads payloads;
  payloads.detector = std::make_unique<detray::io::detector_payload>();
  detray::io::detector_payload& detPayload = *payloads.detector;

  payloads.materialGrids = std::make_unique<detray::io::detector_grids_payload<
      detray::io::surface_material_payload, detray::io::material_id>>();

  detray::io::detector_grids_payload<detray::io::surface_material_payload,
                                     detray::io::material_id>& materialGrids =
      *payloads.materialGrids;

  payloads.surfaceGrids = std::make_unique<
      detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>>();

  detray::io::detector_grids_payload<std::size_t, detray::io::accel_id>&
      surfaceGrids = *payloads.surfaceGrids;

  std::unordered_map<const TrackingVolume*, std::size_t> volumeIds;

  auto lookup = [&volumeIds](const TrackingVolume* v) {
    return volumeIds.at(v);
  };

  std::unordered_map<const TrackingVolume*,
                     std::unordered_map<const Surface*, std::size_t>>
      volumeSurfaceIndices;

  geometry.apply([&](const TrackingVolume& volume) {
    auto& volPayload =
        detPayload.volumes.emplace_back(convertVolume(gctx, volume));
    volPayload.index.link = detPayload.volumes.size() - 1;
    volumeIds[&volume] = volPayload.index.link;

    ACTS_DEBUG("Volume " << volume.volumeName() << " has index "
                         << volPayload.index.link);

    auto& surfaceIndices = volumeSurfaceIndices[&volume];

    for (auto& surface : volume.surfaces()) {
      auto& srfPayload =
          volPayload.surfaces.emplace_back(convertSurface(gctx, surface));
      srfPayload.index_in_coll = volPayload.surfaces.size() - 1;
      srfPayload.masks.at(0).volume_link.link = volPayload.index.link;
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

    ACTS_DEBUG("Volume " << volume.volumeName() << " (detray idx: "
                         << volPayload.index.link << ") has "
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
                              << ", type: " << surface->type()
                              << ") has detray index " << idx);
    }

    // Portals have produced surfaces and are added in volume payload, handle
    // material now

    auto [grids, homogeneous] =
        convertMaterial(volume, surfaceIndices, volPayload);

    ACTS_DEBUG("Volume " << volume.volumeName()
                         << " (detray idx: " << volPayload.index.link
                         << ") has " << homogeneous.surface_mat.size()
                         << " material slabs");

    if (!homogeneous.surface_mat.empty()) {
      // Only add if it's not empty (it might be)
      // NOTE: Currently, it'll always be populated by at least the homogeneous
      // NOTE: Volume association is internal to
      // `detray::io::material_volume_payload`
      if (!payloads.homogeneousMaterial) {
        payloads.homogeneousMaterial = std::make_unique<
            detray::io::detector_homogeneous_material_payload>();
      }
      payloads.homogeneousMaterial->volumes.emplace_back(
          std::move(homogeneous));
    }

    ACTS_DEBUG("Volume " << volume.volumeName()
                         << " (detray idx: " << volPayload.index.link
                         << ") has " << grids.size() << " material grids");
    if (!grids.empty()) {
      // Only add if we have grids
      // NOTE: Volume association is EXTERNAL, i.e. we need to fill a map keyed
      // by the volume index
      materialGrids.grids[volPayload.index.link] = std::move(grids);
    }

    // Look for navigation policies that we need to convert!
    const auto* navPolicy = volume.navigationPolicy();
    if (navPolicy != nullptr) {
      // Create surface lookup function for this volume
      auto surfaceLookupFn =
          [&surfaceIndices](const Surface* surface) -> std::size_t {
        auto it = surfaceIndices.find(surface);
        if (it == surfaceIndices.end()) {
          throw std::runtime_error("Surface not found in surface indices map");
        }
        return it->second;
      };

      std::optional<DetraySurfaceGrid> detrayGrid = std::nullopt;

      navPolicy->visit([&](const INavigationPolicy& policy) {
        auto grid = m_cfg.convertNavigationPolicy(policy, gctx, surfaceLookupFn,
                                                  logger());
        if (!grid.has_value()) {
          // Policies without an explicit detray conversion (see
          // NOOP_CONVERTER_IMPL) are not a conflict: a volume may legitimately
          // combine e.g. a SurfaceArrayNavigationPolicy with a
          // TryAllNavigationPolicy for passives and portals.
          return;
        }

        if (detrayGrid.has_value()) {
          ACTS_ERROR("Volume "
                     << volume.volumeName()
                     << " has more than one detray-convertible navigation "
                        "policy. This cannot currently be handled.");
          throw std::runtime_error{
              "Multiple detray-compatible navigation policies"};
        }

        detrayGrid = std::move(grid);
      });

      if (detrayGrid.has_value()) {
        ACTS_DEBUG("Volume " << volume.volumeName()
                             << " (detray idx: " << volPayload.index.link
                             << ") has navigation policy which produced "
                             << detrayGrid->bins.size() << " populated bins");

        detrayGrid->owner_link.link = volPayload.index.link;

        // Add the surface grid to the payload
        surfaceGrids.grids[volPayload.index.link].push_back(*detrayGrid);
      }
      // per volume, we have a VECTOR of grids: what are they? are they always
      // tied to a surface? which one?
    }
  });

  // HACK: Beampipe MUST have index 0
  std::size_t beampipeIdx = volumeIds.at(m_cfg.beampipeVolume);
  ACTS_DEBUG("Beampipe volume (" << m_cfg.beampipeVolume->volumeName()
                                 << ") index: " << beampipeIdx);
  ACTS_DEBUG("Volume at index 0 is " << detPayload.volumes.at(0).name);

  // Swap beampipe and world volumes self-index
  std::swap(detPayload.volumes.at(0).index.link,
            detPayload.volumes.at(beampipeIdx).index.link);

  // Swap beampipe and world volumes location in vector
  std::swap(detPayload.volumes.at(0), detPayload.volumes.at(beampipeIdx));

  // Adjust volume indices in surfaces after swapping
  for (auto& vol : detPayload.volumes) {
    for (auto& srf : vol.surfaces) {
      // Portals can carry one mask per neighbour volume
      for (auto& mask : srf.masks) {
        if (mask.volume_link.link == beampipeIdx) {
          mask.volume_link.link = 0;
        } else if (mask.volume_link.link == 0) {
          mask.volume_link.link = beampipeIdx;
        }
      }
    }
  }

  if (payloads.homogeneousMaterial) {
    ACTS_DEBUG("Adjusting homogeneous material entries after swapping");
    auto& dthmPayload = *payloads.homogeneousMaterial;

    // Possibly swap homogeneous material entries in vector if they both exist
    auto find = [](std::size_t id) {
      return [id](const auto& vol) { return vol.volume_link.link == id; };
    };

    auto beampipeIt =
        std::ranges::find_if(dthmPayload.volumes, find(beampipeIdx));
    auto worldIt = std::ranges::find_if(dthmPayload.volumes, find(0));

    if (beampipeIt != dthmPayload.volumes.end() &&
        worldIt != dthmPayload.volumes.end()) {
      // BOTH world and beampipe have homogoenous material: swap them
      ACTS_DEBUG("Swapping beampipe and world homogoenous material entries");
      std::swap(*beampipeIt, *worldIt);
    }

    // Retarget the entries, regardless of whether there is an entry for only
    // one of them
    for (auto& mat : dthmPayload.volumes) {
      if (mat.volume_link.link == beampipeIdx) {
        ACTS_DEBUG("Reassigning beampipe homogoenous material to index 0");
        mat.volume_link.link = 0;
      } else if (mat.volume_link.link == 0) {
        ACTS_DEBUG("Reassigning world homogoenous material to beampipe index "
                   << beampipeIdx);
        mat.volume_link.link = beampipeIdx;
      }
    }
  } else {
    ACTS_DEBUG("No homogeneous material payload to adjust after swapping");
  }

  {
    // Adjust material grids after swapping
    auto beampipeGridIt = materialGrids.grids.find(beampipeIdx);
    auto worldGridIt = materialGrids.grids.find(0);

    if (beampipeGridIt != materialGrids.grids.end() &&
        worldGridIt != materialGrids.grids.end()) {
      // BOTH world and beampipe have grid specifiers: swap them
      ACTS_DEBUG("Swapping beampipe and world material grid specifiers");
      std::swap(beampipeGridIt->second, worldGridIt->second);
    } else if (beampipeGridIt != materialGrids.grids.end()) {
      // ONLY beampipe has grid specifier: move it to world
      ACTS_DEBUG("Moving beampipe material grid specifier to world");
      materialGrids.grids[0] = std::move(beampipeGridIt->second);
      materialGrids.grids.erase(beampipeGridIt);
    } else if (worldGridIt != materialGrids.grids.end()) {
      // ONLY world has grid specifier: move it to beampipe
      ACTS_DEBUG("Moving world material grid specifier to beampipe");
      materialGrids.grids[beampipeIdx] = std::move(worldGridIt->second);
      materialGrids.grids.erase(worldGridIt);
    }
  }

  {
    // Adjust surface grids after swapping
    // @NOTE: The beampipe should **generally** not have a surface grid, but
    //        let's be safe and swap them regardless

    auto beampipeGridIt = surfaceGrids.grids.find(beampipeIdx);
    auto worldGridIt = surfaceGrids.grids.find(0);

    if (beampipeGridIt != surfaceGrids.grids.end() &&
        worldGridIt != surfaceGrids.grids.end()) {
      // BOTH world and beampipe have grid specifiers: swap them
      ACTS_DEBUG("Swapping beampipe and world surface grid specifiers");
      std::swap(beampipeGridIt->second, worldGridIt->second);
    } else if (beampipeGridIt != surfaceGrids.grids.end()) {
      // ONLY beampipe has grid specifier: move it to world
      ACTS_DEBUG("Moving beampipe surface grid specifier to world");
      surfaceGrids.grids[0] = std::move(beampipeGridIt->second);
      surfaceGrids.grids.erase(beampipeGridIt);
    } else if (worldGridIt != surfaceGrids.grids.end()) {
      // ONLY world has grid specifier: move it to beampipe
      ACTS_DEBUG("Moving world surface grid specifier to beampipe");
      surfaceGrids.grids[beampipeIdx] = std::move(worldGridIt->second);
      surfaceGrids.grids.erase(worldGridIt);
    }
  }

  // This needs to happen after swapping so that the indices are correct
  payloads.names = {{0, "Detector"}};
  for (const auto& volume : detPayload.volumes) {
    payloads.names.emplace(volume.index.link + 1, volume.name);
  }

  ACTS_DEBUG("Collected " << detPayload.volumes.size() << " volumes");

  return payloads;
}

}  // namespace ActsPlugins
