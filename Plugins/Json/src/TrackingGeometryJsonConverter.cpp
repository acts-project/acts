// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

#include "Acts/Geometry/CompositePortalLink.hpp"
#include "Acts/Geometry/ConeVolumeBounds.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/CutoutCylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GenericCuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrapezoidVolumeBounds.hpp"
#include "Acts/Geometry/TrivialPortalLink.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiNavigationPolicy.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/SurfacePlacementBase.hpp"
#include "Acts/Utilities/AnyGridView.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/AlgebraJsonConverter.hpp"
#include "ActsPlugins/Json/GeometryIdentifierJsonConverter.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <stdexcept>
#include <unordered_set>

template <typename object_t, const char* kContext>
struct Acts::TrackingGeometryJsonConverter::PointerToIdLookup {
  // Insert a new object to ID mapping.
  //
  // @param object is the source object pointer key
  // @param objectId is the serialized ID to assign
  //
  // @return true if insertion happened, false if the object was already
  //         present
  bool emplace(const object_t& object, std::size_t objectId) {
    return m_objectIds.emplace(&object, objectId).second;
  }

  // Resolve a serialized object ID from an object reference.
  //
  // @param object is the source object key
  //
  // @return associated serialized object ID
  //
  // @throw std::invalid_argument if the object is not in the lookup
  std::size_t at(const object_t& object) const {
    auto it = m_objectIds.find(&object);
    if (it == m_objectIds.end()) {
      throw std::invalid_argument("Pointer-to-ID lookup failed for " +
                                  std::string{kContext} +
                                  ": object is outside serialized hierarchy");
    }
    return it->second;
  }

 private:
  // Container mapping object pointers to their respective id
  std::unordered_map<const object_t*, std::size_t> m_objectIds;
};

template <typename object_t, typename pointer_t, const char* kContext>
struct Acts::TrackingGeometryJsonConverter::IdToPointerLikeLookup {
  /// Insert a new serialized ID to object mapping.
  ///
  /// @param objectId is the serialized ID key
  /// @param object is the target pointer-like object
  ///
  /// @return true if insertion happened, false if the ID was already present
  bool emplace(std::size_t objectId, pointer_t object) {
    return m_objects.emplace(objectId, std::move(object)).second;
  }

  /// Try to find a mapped pointer-like object by serialized ID.
  ///
  /// @param objectId is the serialized ID key
  ///
  /// @return mapped pointer-like object, or null-equivalent if not found
  pointer_t find(std::size_t objectId) const {
    auto it = m_objects.find(objectId);
    return it == m_objects.end() ? pointer_t{} : it->second;
  }

  /// Resolve a mapped pointer-like object by serialized ID.
  ///
  /// @param objectId is the serialized ID key
  ///
  /// @return mapped pointer-like object reference
  ///
  /// @throw std::invalid_argument if the ID is not mapped
  const pointer_t& at(std::size_t objectId) const {
    auto it = m_objects.find(objectId);
    if (it == m_objects.end()) {
      throw std::invalid_argument("ID-to-pointer lookup failed for " +
                                  std::string{kContext} +
                                  ": unknown serialized object ID");
    }
    return it->second;
  }

 private:
  /// Container mapping object ids to their respective pointers
  std::unordered_map<std::size_t, pointer_t> m_objects;
};

namespace {

// Internal notes on the conversion flow used in this translation unit:
// - Encoder side:
//   - normalize enum-like values (axis/shape metadata) to stable JSON strings
//   - encode bounds and portal-link polymorphic payloads using kind tags
//   - emit volumes in deterministic depth-first order and connect them by IDs
//   - emit unique portals in a top-level table and refer to them from volumes
// - Decoder side:
//   - read and validate schema metadata
//   - construct all volumes first, then attach children and portals
//   - rebuild polymorphic portal links from kind-tagged payloads

constexpr const char* kHeaderKey = "acts-geometry-volumes";
constexpr const char* kVersionKey = "format-version";
constexpr const char* kScopeKey = "scope";
constexpr const char* kScopeValue = "volumes-bounds-portals";
constexpr int kFormatVersion = 1;

constexpr const char* kPortalsKey = "portals";
constexpr const char* kSurfacesKey = "surfaces";
constexpr const char* kVolumesKey = "volumes";

constexpr const char* kRootVolumeIdKey = "root_volume_id";
constexpr const char* kVolumeIdKey = "volume_id";
constexpr const char* kPortalIdKey = "portal_id";
constexpr const char* kSurfaceIdKey = "surface_id";

constexpr const char* kNameKey = "name";
constexpr const char* kGeometryIdKey = "geometry_id";
constexpr const char* kTransformKey = "transform";
constexpr const char* kBoundsKey = "bounds";
constexpr const char* kChildrenKey = "children";

constexpr const char* kPortalIdsKey = "portal_ids";
constexpr const char* kAlongNormalKey = "along_normal";
constexpr const char* kOppositeNormalKey = "opposite_normal";

constexpr const char* kNavigationPolicyKey = "navigation_policy";

constexpr const char* kKindKey = "kind";
constexpr const char* kValuesKey = "values";
constexpr const char* kTargetVolumeIdKey = "target_volume_id";
constexpr const char* kDirectionKey = "direction";
constexpr const char* kAxesKey = "axes";
constexpr const char* kBinsKey = "bins";
constexpr const char* kLocalKey = "local";
constexpr const char* kArtifactLinksKey = "artifact_links";

// -------------------------------------------------------------------
// Kind getters

template <typename bounds_t>
std::string getPortalLinkKind() {
  if (std::is_same_v<bounds_t, Acts::TrivialPortalLink>) {
    return "Trivial";
  } else if (std::is_same_v<bounds_t, Acts::CompositePortalLink>) {
    return "Composite";
  } else if (std::is_same_v<bounds_t, Acts::GridPortalLink>) {
    return "Grid";
  } else {
    throw std::invalid_argument("Unknown portal link kind");
  }
}

template <typename bounds_t>
std::string getVolumeBoundsKind() {
  if (std::is_same_v<bounds_t, Acts::ConeVolumeBounds>) {
    return "Cone";
  } else if (std::is_same_v<bounds_t, Acts::CuboidVolumeBounds>) {
    return "Cuboid";
  } else if (std::is_same_v<bounds_t, Acts::CutoutCylinderVolumeBounds>) {
    return "CutoutCylinder";
  } else if (std::is_same_v<bounds_t, Acts::CylinderVolumeBounds>) {
    return "Cylinder";
  } else if (std::is_same_v<bounds_t, Acts::GenericCuboidVolumeBounds>) {
    return "GenericCuboid";
  } else if (std::is_same_v<bounds_t, Acts::TrapezoidVolumeBounds>) {
    return "Trapezoid";
  } else {
    throw std::invalid_argument("Unknown volume bounds kind");
  }
}

template <typename bounds_t>
std::string getNavigationPolicyKind() {
  if (std::is_same_v<bounds_t, Acts::TryAllNavigationPolicy>) {
    return "TryAll";
  } else if (std::is_same_v<bounds_t, Acts::SurfaceArrayNavigationPolicy>) {
    return "SurfaceArray";
  } else if (std::is_same_v<bounds_t, Acts::MultiNavigationPolicy>) {
    return "MultiNavigation";
  } else {
    throw std::invalid_argument("Unknown portal link kind");
  }
}

// -------------------------------------------------------------------
// Volume bounds encoder/decoder

template <typename bounds_t>
std::unique_ptr<Acts::VolumeBounds> decodeVolumeBoundsT(
    const nlohmann::json& jBounds) {
  constexpr std::size_t kValues = bounds_t::BoundValues::eSize;
  const auto values = jBounds.at(kValuesKey).get<std::vector<double>>();
  if (values.size() != kValues) {
    throw std::invalid_argument("Invalid number of values for volume bounds");
  }
  std::array<double, kValues> boundValues{};
  std::copy_n(values.begin(), kValues, boundValues.begin());
  return std::make_unique<bounds_t>(boundValues);
}

template <typename bounds_t>
nlohmann::json encodeVolumeBoundsT(const bounds_t& bounds) {
  nlohmann::json jBounds;
  jBounds["kind"] = getVolumeBoundsKind<bounds_t>();
  jBounds[kValuesKey] = bounds.values();
  return jBounds;
}

// -------------------------------------------------------------------
// Navigation policy encoder/decoder

std::unique_ptr<Acts::INavigationPolicy> decodeMultiNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  std::vector<std::unique_ptr<Acts::INavigationPolicy>> children;
  for (const auto& child : encoded.at(kChildrenKey)) {
    children.push_back(converter.navigationPolicyFromJson(
        gctx, child[kNavigationPolicyKey], volume, logger));
  }
  return std::make_unique<Acts::MultiNavigationPolicy>(std::move(children));
}

std::unique_ptr<Acts::INavigationPolicy> decodeTryAllNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  Acts::TryAllNavigationPolicy::Config cfg;
  cfg.passives = encoded.at("passives").get<bool>();
  cfg.sensitives = encoded.at("sensitives").get<bool>();
  cfg.portals = encoded.at("portals").get<bool>();

  return std::make_unique<Acts::TryAllNavigationPolicy>(gctx, volume, logger,
                                                        cfg);
}

std::unique_ptr<Acts::INavigationPolicy> decodeSurfaceArrayNavigationPolicy(
    const nlohmann::json& encoded, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) {
  Acts::SurfaceArrayNavigationPolicy::Config cfg;
  cfg.layerType = encoded.at("layerType")
                      .get<Acts::SurfaceArrayNavigationPolicy::LayerType>();
  cfg.bins = {encoded.at("bins0").get<std::size_t>(),
              encoded.at("bins1").get<std::size_t>()};

  return std::make_unique<Acts::SurfaceArrayNavigationPolicy>(gctx, volume,
                                                              logger, cfg);
}

nlohmann::json encodeMultiNavigationPolicy(
    const Acts::MultiNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& converter) {
  nlohmann::json jPolicy;
  jPolicy[kKindKey] = getNavigationPolicyKind<Acts::MultiNavigationPolicy>();
  for (const auto& pol : policy.policies()) {
    nlohmann::json jPol;
    jPol["navigation_policy"] = converter.navigationPolicyToJson(*pol);
    jPolicy[kChildrenKey].push_back(jPol);
  }
  return jPolicy;
}

nlohmann::json encodeTryAllNavigationPolicy(
    const Acts::TryAllNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& /*converter*/) {
  const auto& cfg = policy.config();

  nlohmann::json jPolicy;
  jPolicy[kKindKey] = getNavigationPolicyKind<Acts::TryAllNavigationPolicy>();
  jPolicy["portals"] = cfg.portals;
  jPolicy["sensitives"] = cfg.sensitives;
  jPolicy["passives"] = cfg.passives;
  return jPolicy;
}

nlohmann::json encodeSurfaceArrayNavigationPolicy(
    const Acts::SurfaceArrayNavigationPolicy& policy,
    const Acts::TrackingGeometryJsonConverter& /*converter*/) {
  const auto& cfg = policy.config();

  nlohmann::json jPolicy;
  jPolicy[kKindKey] =
      getNavigationPolicyKind<Acts::SurfaceArrayNavigationPolicy>();
  jPolicy["layerType"] = cfg.layerType;
  jPolicy["bins0"] = cfg.bins.first;
  jPolicy["bins1"] = cfg.bins.second;
  return jPolicy;
}

// -------------------------------------------------------------------
// Portal link encoder/decoder

std::shared_ptr<Acts::RegularSurface> regularSurfaceFromJson(
    const nlohmann::json& jSurface) {
  auto surface = Acts::SurfaceJsonConverter::fromJson(jSurface);
  auto regular = std::dynamic_pointer_cast<Acts::RegularSurface>(surface);
  if (regular == nullptr) {
    throw std::invalid_argument("Portal link surface is not a RegularSurface");
  }
  return regular;
}

std::unique_ptr<Acts::GridPortalLink> makeGridPortalLink(
    const std::shared_ptr<Acts::RegularSurface>& surface,
    Acts::AxisDirection direction, const Acts::IAxis& axis0,
    const Acts::IAxis* axis1) {
  std::unique_ptr<Acts::GridPortalLink> grid;

  if (axis1 == nullptr) {
    axis0.visit([&](const auto& a0) {
      using axis_t = std::remove_cvref_t<decltype(a0)>;
      axis_t axisCopy = a0;
      grid =
          Acts::GridPortalLink::make(surface, direction, std::move(axisCopy));
    });
  } else {
    axis0.visit([&](const auto& a0) {
      using axis0_t = std::remove_cvref_t<decltype(a0)>;
      axis0_t axis0Copy = a0;
      axis1->visit([&](const auto& a1) {
        using axis1_t = std::remove_cvref_t<decltype(a1)>;
        axis1_t axis1Copy = a1;
        grid = Acts::GridPortalLink::make(surface, std::move(axis0Copy),
                                          std::move(axis1Copy));
      });
    });
  }

  if (grid == nullptr) {
    throw std::invalid_argument("Could not construct GridPortalLink from axes");
  }

  if (grid->direction() != direction) {
    throw std::invalid_argument(
        "Decoded grid direction does not match payload");
  }

  return grid;
}

nlohmann::json encodeTrivialPortalLink(
    const Acts::TrivialPortalLink& link, const Acts::GeometryContext& /*gctx*/,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = getPortalLinkKind<Acts::TrivialPortalLink>();
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kTargetVolumeIdKey] = volumeIds.at(link.volume());
  return jLink;
}

nlohmann::json encodeCompositePortalLink(
    const Acts::CompositePortalLink& link, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = getPortalLinkKind<Acts::CompositePortalLink>();
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kDirectionKey] = link.direction();
  jLink[kChildrenKey] = nlohmann::json::array();

  for (const auto& child : link.links()) {
    jLink[kChildrenKey].push_back(
        converter.portalLinkToJson(gctx, child, surfaceIds, volumeIds));
  }
  return jLink;
}

nlohmann::json encodeGridPortalLink(
    const Acts::GridPortalLink& link, const Acts::GeometryContext& gctx,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    const Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  nlohmann::json jLink;
  jLink[kKindKey] = getPortalLinkKind<Acts::GridPortalLink>();
  jLink[kDirectionKey] = link.direction();
  jLink[kSurfaceIdKey] = surfaceIds.at(link.surface());
  jLink[kAxesKey] = nlohmann::json::array();

  for (const auto* axis : link.grid().axes()) {
    jLink[kAxesKey].push_back(Acts::AxisJsonConverter::toJson(*axis));
  }

  Acts::AnyGridConstView<const Acts::TrackingVolume*> view(link.grid());
  const auto nBins = view.numLocalBins();
  const auto dim = view.dimensions();

  jLink[kBinsKey] = nlohmann::json::array();
  if (dim == 1u) {
    for (std::size_t i0 = 0u; i0 <= nBins.at(0) + 1u; ++i0) {
      nlohmann::json jBin;
      jBin[kLocalKey] = std::vector<std::size_t>{i0};
      if (const auto* target = view.atLocalBins({i0}); target == nullptr) {
        jBin[kTargetVolumeIdKey] = nullptr;
      } else {
        jBin[kTargetVolumeIdKey] = volumeIds.at(*target);
      }
      jLink[kBinsKey].push_back(std::move(jBin));
    }
  } else if (dim == 2u) {
    for (std::size_t i0 = 0u; i0 <= nBins.at(0) + 1u; ++i0) {
      for (std::size_t i1 = 0u; i1 <= nBins.at(1) + 1u; ++i1) {
        nlohmann::json jBin;
        jBin[kLocalKey] = std::vector<std::size_t>{i0, i1};
        if (const auto* target = view.atLocalBins({i0, i1});
            target == nullptr) {
          jBin[kTargetVolumeIdKey] = nullptr;
        } else {
          jBin[kTargetVolumeIdKey] = volumeIds.at(*target);
        }
        jLink[kBinsKey].push_back(std::move(jBin));
      }
    }
  } else {
    throw std::invalid_argument("Unsupported GridPortalLink dimensionality");
  }

  jLink[kArtifactLinksKey] = nlohmann::json::array();
  for (const auto& artifact : link.artifactPortalLinks()) {
    jLink[kArtifactLinksKey].push_back(
        converter.portalLinkToJson(gctx, artifact, surfaceIds, volumeIds));
  }

  return jLink;
}

std::unique_ptr<Acts::PortalLinkBase> decodeTrivialPortalLink(
    const nlohmann::json& encoded,
    const Acts::TrackingGeometryJsonConverter& /*converter*/,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  const auto linkSurfaceId = encoded.at(kSurfaceIdKey).get<std::size_t>();
  const auto targetVolumeId = encoded.at(kTargetVolumeIdKey).get<std::size_t>();
  return std::make_unique<Acts::TrivialPortalLink>(surfaces.at(linkSurfaceId),
                                                   *volumes.at(targetVolumeId));
}

std::unique_ptr<Acts::PortalLinkBase> decodeCompositePortalLink(
    const nlohmann::json& encoded,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  const auto direction = encoded.at(kDirectionKey).get<Acts::AxisDirection>();
  std::vector<std::unique_ptr<Acts::PortalLinkBase>> children;
  for (const auto& child : encoded.at(kChildrenKey)) {
    children.push_back(converter.portalLinkFromJson(child, surfaces, volumes));
  }
  return std::make_unique<Acts::CompositePortalLink>(std::move(children),
                                                     direction);
}

std::unique_ptr<Acts::PortalLinkBase> decodeGridPortalLink(
    const nlohmann::json& encoded,
    const Acts::TrackingGeometryJsonConverter& converter,
    const Acts::TrackingGeometryJsonConverter::SurfacePointerLookup& surfaces,
    const Acts::TrackingGeometryJsonConverter::VolumePointerLookup& volumes) {
  auto linkSurfaceId = encoded.at(kSurfaceIdKey).get<std::size_t>();
  const auto direction = encoded.at(kDirectionKey).get<Acts::AxisDirection>();

  std::vector<std::unique_ptr<Acts::IAxis>> axes;
  for (const auto& jAxis : encoded.at(kAxesKey)) {
    axes.push_back(Acts::AxisJsonConverter::fromJson(jAxis));
  }
  if (axes.empty() || axes.size() > 2u) {
    throw std::invalid_argument("GridPortalLink requires 1 or 2 axes");
  }

  auto grid =
      makeGridPortalLink(surfaces.at(linkSurfaceId), direction, *axes.at(0),
                         axes.size() == 2u ? axes.at(1).get() : nullptr);

  Acts::AnyGridView<const Acts::TrackingVolume*> view(grid->grid());
  const auto dim = view.dimensions();
  for (const auto& jBin : encoded.at(kBinsKey)) {
    auto local = jBin.at(kLocalKey).get<std::vector<std::size_t>>();
    if (local.size() != dim) {
      throw std::invalid_argument("Grid bin dimensionality mismatch");
    }
    Acts::IGrid::AnyIndexType localIndices(local.begin(), local.end());

    const Acts::TrackingVolume* target = nullptr;
    if (!jBin.at(kTargetVolumeIdKey).is_null()) {
      const auto targetVolumeId =
          jBin.at(kTargetVolumeIdKey).get<std::size_t>();
      target = volumes.at(targetVolumeId);
    }
    view.atLocalBins(localIndices) = target;
  }

  std::vector<Acts::TrivialPortalLink> artifacts;
  if (encoded.contains(kArtifactLinksKey)) {
    for (const auto& jArtifact : encoded.at(kArtifactLinksKey)) {
      auto decodedArtifact =
          converter.portalLinkFromJson(jArtifact, surfaces, volumes);
      auto* trivial =
          dynamic_cast<Acts::TrivialPortalLink*>(decodedArtifact.get());
      if (trivial == nullptr) {
        throw std::invalid_argument(
            "GridPortalLink artifact link is not trivial");
      }
      artifacts.push_back(std::move(*trivial));
    }
  }
  grid->setArtifactPortalLinks(std::move(artifacts));

  return grid;
}

// -------------------------------------------------------------------
// Records for temporary storage

struct SurfaceRecord {
  std::size_t surfaceId = 0u;
  nlohmann::json payload;
};

struct PortalRecord {
  std::size_t portalId = 0u;
  nlohmann::json payload;
};

struct VolumeRecord {
  std::size_t volumeId = 0u;
  std::string name;
  Acts::GeometryIdentifier::Value geometryId = 0u;
  Acts::Transform3 transform{Acts::Transform3::Identity()};
  nlohmann::json bounds;
  std::vector<std::size_t> children;
  std::vector<std::size_t> portalIds;
  std::vector<std::size_t> surfaceIds;
  nlohmann::json navigationPolicy;
};

// -------------------------------------------------------------------
// Utilities

void verifySchemaHeader(const nlohmann::json& encoded) {
  if (!encoded.contains(kHeaderKey)) {
    throw std::invalid_argument("Missing geometry JSON header");
  }
  const auto& header = encoded.at(kHeaderKey);
  if (header.at(kVersionKey).get<int>() != kFormatVersion) {
    throw std::invalid_argument("Unsupported geometry JSON format version");
  }
  if (header.at(kScopeKey).get<std::string>() != kScopeValue) {
    throw std::invalid_argument("Unexpected geometry JSON scope");
  }
}

void ensureIdentifiers(Acts::TrackingVolume& volume,
                       Acts::GeometryIdentifier::Value& nextVolumeId) {
  Acts::GeometryIdentifier volumeId = volume.geometryId();
  if (volumeId == Acts::GeometryIdentifier{}) {
    volumeId = Acts::GeometryIdentifier{}.withVolume(nextVolumeId++);
    volume.assignGeometryId(volumeId);
  }

  for (const auto [ib, boundary] : Acts::enumerate(volume.boundarySurfaces())) {
    // Gen1 api, ignore Sonar complaints
    auto& mutableBoundarySurface =
        const_cast<Acts::RegularSurface&>(boundary->surfaceRepresentation());
    if (mutableBoundarySurface.geometryId() == Acts::GeometryIdentifier{}) {
      mutableBoundarySurface.assignGeometryId(
          Acts::GeometryIdentifier(volumeId).withBoundary(ib + 1u));
    }
  }

  std::size_t portalIndex = 0u;
  for (auto& portal : volume.portals()) {
    if (auto& mutablePortalSurface = portal.surface();
        mutablePortalSurface.geometryId() == Acts::GeometryIdentifier{}) {
      mutablePortalSurface.assignGeometryId(
          Acts::GeometryIdentifier(volumeId).withExtra(portalIndex + 1u));
    }
    ++portalIndex;
  }

  for (auto& child : volume.volumes()) {
    ensureIdentifiers(child, nextVolumeId);
  }
}

void collectGeometry(
    const Acts::TrackingVolume& volume,
    std::vector<const Acts::Surface*>& orderedSurfaces,
    std::vector<const Acts::Portal*>& orderedPortals,
    std::vector<const Acts::TrackingVolume*>& orderedVolumes,
    Acts::TrackingGeometryJsonConverter::SurfaceIdLookup& surfaceIds,
    Acts::TrackingGeometryJsonConverter::PortalIdLookup& portalIds,
    Acts::TrackingGeometryJsonConverter::VolumeIdLookup& volumeIds) {
  auto insertSurface = [&](const auto& surf) {
    if (surfaceIds.emplace(surf, orderedSurfaces.size())) {
      orderedSurfaces.push_back(&surf);
    }
  };
  auto insertPortal = [&](const auto& port) {
    if (portalIds.emplace(port, orderedPortals.size())) {
      orderedPortals.push_back(&port);
    }
  };
  auto insertVolume = [&](const auto& vol) {
    if (volumeIds.emplace(vol, orderedVolumes.size())) {
      orderedVolumes.push_back(&vol);
    }
  };

  auto insertPortalLink = [&](const auto& link, auto&& self) {
    insertSurface(link->surface());

    if (const auto* trivial =
            dynamic_cast<const Acts::TrivialPortalLink*>(link);
        trivial != nullptr) {
      return;
    }

    if (const auto* composite =
            dynamic_cast<const Acts::CompositePortalLink*>(link);
        composite != nullptr) {
      const auto& children = composite->links();
      for (const auto& child : children) {
        self(&child, self);
      }
      return;
    }

    const auto* grid = dynamic_cast<const Acts::GridPortalLink*>(link);
    if (grid != nullptr) {
      const auto& children = grid->artifactPortalLinks();
      for (const auto& child : children) {
        self(&child, self);
      }
      return;
    }
  };

  insertVolume(volume);

  for (const auto& portal : volume.portals()) {
    insertPortal(portal);

    insertSurface(portal.surface());

    if (const auto* along = portal.getLink(Acts::Direction::AlongNormal());
        along != nullptr) {
      insertPortalLink(along, insertPortalLink);
    }

    if (const auto* opposite =
            portal.getLink(Acts::Direction::OppositeNormal());
        opposite != nullptr) {
      insertPortalLink(opposite, insertPortalLink);
    }
  }
  for (const auto& surf : volume.surfaces()) {
    insertSurface(surf);
  }

  for (const auto& child : volume.volumes()) {
    collectGeometry(child, orderedSurfaces, orderedPortals, orderedVolumes,
                    surfaceIds, portalIds, volumeIds);
  }
}

}  // namespace

Acts::TrackingGeometryJsonConverter::Config
Acts::TrackingGeometryJsonConverter::Config::defaultConfig() {
  Config cfg;

  cfg.encodeVolumeBounds.registerFunction(encodeVolumeBoundsT<ConeVolumeBounds>)
      .registerFunction(encodeVolumeBoundsT<CuboidVolumeBounds>)
      .registerFunction(encodeVolumeBoundsT<CutoutCylinderVolumeBounds>)
      .registerFunction(encodeVolumeBoundsT<CylinderVolumeBounds>)
      .registerFunction(encodeVolumeBoundsT<GenericCuboidVolumeBounds>)
      .registerFunction(encodeVolumeBoundsT<TrapezoidVolumeBounds>);

  cfg.encodeNavigationPolicy.registerFunction(encodeTryAllNavigationPolicy)
      .registerFunction(encodeSurfaceArrayNavigationPolicy)
      .registerFunction(encodeMultiNavigationPolicy);

  cfg.encodePortalLink.registerFunction(encodeTrivialPortalLink)
      .registerFunction(encodeCompositePortalLink)
      .registerFunction(encodeGridPortalLink);

  cfg.decodePortalLink
      .registerKind(getPortalLinkKind<TrivialPortalLink>(),
                    decodeTrivialPortalLink)
      .registerKind(getPortalLinkKind<CompositePortalLink>(),
                    decodeCompositePortalLink)
      .registerKind(getPortalLinkKind<GridPortalLink>(), decodeGridPortalLink);

  cfg.decodeVolumeBounds
      .registerKind(getVolumeBoundsKind<ConeVolumeBounds>(),
                    decodeVolumeBoundsT<ConeVolumeBounds>)
      .registerKind(getVolumeBoundsKind<CuboidVolumeBounds>(),
                    decodeVolumeBoundsT<CuboidVolumeBounds>)
      .registerKind(getVolumeBoundsKind<CutoutCylinderVolumeBounds>(),
                    decodeVolumeBoundsT<CutoutCylinderVolumeBounds>)
      .registerKind(getVolumeBoundsKind<CylinderVolumeBounds>(),
                    decodeVolumeBoundsT<CylinderVolumeBounds>)
      .registerKind(getVolumeBoundsKind<GenericCuboidVolumeBounds>(),
                    decodeVolumeBoundsT<GenericCuboidVolumeBounds>)
      .registerKind(getVolumeBoundsKind<TrapezoidVolumeBounds>(),
                    decodeVolumeBoundsT<TrapezoidVolumeBounds>);

  cfg.decodeNavigationPolicy
      .registerKind(getNavigationPolicyKind<TryAllNavigationPolicy>(),
                    decodeTryAllNavigationPolicy)
      .registerKind(getNavigationPolicyKind<SurfaceArrayNavigationPolicy>(),
                    decodeSurfaceArrayNavigationPolicy)
      .registerKind(getNavigationPolicyKind<MultiNavigationPolicy>(),
                    decodeMultiNavigationPolicy);

  return cfg;
}

Acts::TrackingGeometryJsonConverter::TrackingGeometryJsonConverter(
    Config config)
    : m_cfg(std::move(config)) {}

nlohmann::json Acts::TrackingGeometryJsonConverter::toJson(
    const GeometryContext& gctx, const TrackingGeometry& geometry,
    const Options& options) const {
  if (geometry.geometryVersion() != TrackingGeometry::GeometryVersion::Gen3) {
    throw std::invalid_argument(
        "Tracking geometry serialization is only implemented for Gen3 "
        "geometries");
  }
  return trackingVolumeToJson(gctx, *geometry.highestTrackingVolume(), options);
}

nlohmann::json Acts::TrackingGeometryJsonConverter::portalLinkToJson(
    const GeometryContext& gctx, const PortalLinkBase& link,
    const SurfaceIdLookup& surfaceIds, const VolumeIdLookup& volumeIds) const {
  return m_cfg.encodePortalLink(link, gctx, *this, surfaceIds, volumeIds);
}

std::unique_ptr<Acts::PortalLinkBase>
Acts::TrackingGeometryJsonConverter::portalLinkFromJson(
    const nlohmann::json& encoded, const SurfacePointerLookup& surfaces,
    const VolumePointerLookup& volumes) const {
  return m_cfg.decodePortalLink(encoded, *this, surfaces, volumes);
}

nlohmann::json Acts::TrackingGeometryJsonConverter::navigationPolicyToJson(
    const Acts::INavigationPolicy& policy) const {
  return m_cfg.encodeNavigationPolicy(policy, *this);
}

std::unique_ptr<Acts::INavigationPolicy>
Acts::TrackingGeometryJsonConverter::navigationPolicyFromJson(
    const Acts::GeometryContext& gctx, const nlohmann::json& encoded,
    const Acts::TrackingVolume& volume, const Acts::Logger& logger) const {
  return m_cfg.decodeNavigationPolicy(encoded, gctx, *this, volume, logger);
}

nlohmann::json Acts::TrackingGeometryJsonConverter::trackingVolumeToJson(
    const GeometryContext& gctx, const TrackingVolume& world,
    const Options& /*options*/) const {
  nlohmann::json encoded;
  encoded[kHeaderKey] = nlohmann::json::object();
  encoded[kHeaderKey][kVersionKey] = kFormatVersion;
  encoded[kHeaderKey][kScopeKey] = kScopeValue;

  // Collect object
  std::vector<const Surface*> orderedSurfaces;
  std::vector<const Portal*> orderedPortals;
  std::vector<const TrackingVolume*> orderedVolumes;

  SurfaceIdLookup surfaceIds;
  PortalIdLookup portalIds;
  VolumeIdLookup volumeIds;
  collectGeometry(world, orderedSurfaces, orderedPortals, orderedVolumes,
                  surfaceIds, portalIds, volumeIds);

  encoded[kSurfacesKey] = nlohmann::json::array();
  encoded[kPortalsKey] = nlohmann::json::array();
  encoded[kVolumesKey] = nlohmann::json::array();

  encoded[kRootVolumeIdKey] = volumeIds.at(world);

  // Encode surfaces
  for (const auto* surf : orderedSurfaces) {
    nlohmann::json jSurface = SurfaceJsonConverter::toJson(gctx, *surf);
    jSurface[kSurfaceIdKey] = surfaceIds.at(*surf);
    encoded[kSurfacesKey].push_back(std::move(jSurface));
  }

  // Encode portals
  for (const auto* portal : orderedPortals) {
    nlohmann::json jPortal;
    jPortal[kPortalIdKey] = portalIds.at(*portal);

    if (const auto* along = portal->getLink(Direction::AlongNormal());
        along != nullptr) {
      jPortal[kAlongNormalKey] =
          portalLinkToJson(gctx, *along, surfaceIds, volumeIds);
    } else {
      jPortal[kAlongNormalKey] = nullptr;
    }

    if (const auto* opposite = portal->getLink(Direction::OppositeNormal());
        opposite != nullptr) {
      jPortal[kOppositeNormalKey] =
          portalLinkToJson(gctx, *opposite, surfaceIds, volumeIds);
    } else {
      jPortal[kOppositeNormalKey] = nullptr;
    }

    jPortal[kSurfaceIdKey] = surfaceIds.at(portal->surface());
    encoded[kPortalsKey].push_back(std::move(jPortal));
  }

  // Encode volumes
  for (const auto* volume : orderedVolumes) {
    nlohmann::json jVolume;
    jVolume[kVolumeIdKey] = volumeIds.at(*volume);
    jVolume[kNameKey] = volume->volumeName();
    jVolume[kGeometryIdKey] = nlohmann::json(volume->geometryId());
    jVolume[kTransformKey] =
        Transform3JsonConverter::toJson(volume->localToGlobalTransform(gctx));
    jVolume[kBoundsKey] = m_cfg.encodeVolumeBounds(volume->volumeBounds());

    jVolume[kNavigationPolicyKey] =
        navigationPolicyToJson(*volume->navigationPolicy());

    jVolume[kChildrenKey] = nlohmann::json::array();
    for (const auto& child : volume->volumes()) {
      jVolume[kChildrenKey].push_back(volumeIds.at(child));
    }

    jVolume[kPortalIdsKey] = nlohmann::json::array();
    for (const auto& portal : volume->portals()) {
      jVolume[kPortalIdsKey].push_back(portalIds.at(portal));
    }

    jVolume[kSurfaceIdKey] = nlohmann::json::array();
    for (const auto& surface : volume->surfaces()) {
      jVolume[kSurfaceIdKey].push_back(surfaceIds.at(surface));
    }

    encoded[kVolumesKey].push_back(std::move(jVolume));
  }

  return encoded;
}

std::shared_ptr<Acts::TrackingVolume>
Acts::TrackingGeometryJsonConverter::trackingVolumeFromJson(
    const GeometryContext& gctx, const nlohmann::json& encoded,
    const Options& /*options*/) const {
  verifySchemaHeader(encoded);

  if (!encoded.contains(kVolumesKey) || !encoded.contains(kRootVolumeIdKey) ||
      !encoded.contains(kPortalsKey)) {
    throw std::invalid_argument(
        "Missing volume payload in tracking geometry JSON");
  }

  // Collect surface data
  std::unordered_map<std::size_t, SurfaceRecord> surfaceRecords;
  for (const auto& jSurface : encoded.at(kSurfacesKey)) {
    SurfaceRecord record;
    record.surfaceId = jSurface.at(kSurfaceIdKey).get<std::size_t>();
    record.payload = jSurface;
    const auto [id, inserted] =
        surfaceRecords.try_emplace(record.surfaceId, std::move(record));
    if (!inserted) {
      throw std::invalid_argument("Duplicate serialized surface ID");
    }
  }

  // Collect portal data
  std::unordered_map<std::size_t, PortalRecord> portalRecords;
  for (const auto& jPortal : encoded.at(kPortalsKey)) {
    PortalRecord record;
    record.portalId = jPortal.at(kPortalIdKey).get<std::size_t>();
    record.payload = jPortal;
    const auto [id, inserted] =
        portalRecords.try_emplace(record.portalId, std::move(record));
    if (!inserted) {
      throw std::invalid_argument("Duplicate serialized portal ID");
    }
  }

  // Collect volume data
  std::unordered_map<std::size_t, VolumeRecord> volumeRecords;
  for (const auto& jVolume : encoded.at(kVolumesKey)) {
    VolumeRecord record;
    record.volumeId = jVolume.at(kVolumeIdKey).get<std::size_t>();
    record.name = jVolume.at(kNameKey).get<std::string>();
    record.transform =
        Transform3JsonConverter::fromJson(jVolume.at(kTransformKey));
    record.bounds = jVolume.at(kBoundsKey);
    record.children = jVolume.value(kChildrenKey, std::vector<std::size_t>{});
    record.portalIds = jVolume.value(kPortalIdsKey, std::vector<std::size_t>{});
    record.surfaceIds =
        jVolume.value(kSurfaceIdKey, std::vector<std::size_t>{});
    record.navigationPolicy = jVolume.at(kNavigationPolicyKey);

    if (!jVolume["geometry_id"].is_null()) {
      GeometryIdentifier geoID =
          jVolume["geometry_id"].get<GeometryIdentifier>();
      record.geometryId = geoID.value();
    } else {
      record.geometryId = 0;
    }

    const auto [id, inserted] =
        volumeRecords.try_emplace(record.volumeId, std::move(record));
    if (!inserted) {
      throw std::invalid_argument("Duplicate serialized volume ID");
    }
  }

  // Get root volume id
  const std::size_t rootVolumeId =
      encoded.at(kRootVolumeIdKey).get<std::size_t>();
  if (!volumeRecords.contains(rootVolumeId)) {
    throw std::invalid_argument("Serialized root volume ID does not exist");
  }

  // Collect surface pointers
  SurfacePointerLookup surfacePointers;

  // ---------------------------------------------------
  for (const auto& [surfaceId, record] : surfaceRecords) {
    auto surface = regularSurfaceFromJson(record.payload);
    surfacePointers.emplace(surfaceId, surface);
  }

  // Collect volume pointers
  std::unordered_map<std::size_t, std::unique_ptr<TrackingVolume>>
      volumeStorage;
  VolumePointerLookup volumePointers;

  for (const auto& [volumeId, record] : volumeRecords) {
    auto volumeBounds = m_cfg.decodeVolumeBounds(record.bounds);
    auto volume = std::make_unique<TrackingVolume>(
        record.transform, std::move(volumeBounds), record.name);

    GeometryIdentifier geometryId(record.geometryId);
    if (geometryId == GeometryIdentifier{}) {
      geometryId = GeometryIdentifier{}.withVolume(volumeId + 1u);
    }
    volume->assignGeometryId(geometryId);
    volumePointers.emplace(volumeId, volume.get());
    volumeStorage.emplace(volumeId, std::move(volume));
  }

  std::unordered_set<std::size_t> visiting;
  std::unordered_set<std::size_t> built;

  // Assemble the volume hierarchy
  auto attachChildren = [&](std::size_t volumeId, auto&& self) {
    if (built.contains(volumeId)) {
      return;
    }
    if (!visiting.insert(volumeId).second) {
      throw std::invalid_argument(
          "Cycle detected in serialized volume hierarchy");
    }

    const auto& parent = volumeStorage.at(volumeId);
    if (parent == nullptr) {
      throw std::invalid_argument("Volume was already moved unexpectedly");
    }

    for (std::size_t childId : volumeRecords.at(volumeId).children) {
      self(childId, self);
      auto& child = volumeStorage.at(childId);
      if (child == nullptr) {
        throw std::invalid_argument(
            "Serialized child volume has already been attached");
      }
      parent->addVolume(std::move(child));
    }

    visiting.erase(volumeId);
    built.insert(volumeId);
  };

  attachChildren(rootVolumeId, attachChildren);

  auto decodePortal = [&](const nlohmann::json& jPortal) {
    std::unique_ptr<PortalLinkBase> along = nullptr;
    std::unique_ptr<PortalLinkBase> opposite = nullptr;

    if (!jPortal.at(kAlongNormalKey).is_null()) {
      along = portalLinkFromJson(jPortal.at(kAlongNormalKey), surfacePointers,
                                 volumePointers);
    }
    if (!jPortal.at(kOppositeNormalKey).is_null()) {
      opposite = portalLinkFromJson(jPortal.at(kOppositeNormalKey),
                                    surfacePointers, volumePointers);
    }
    if (along == nullptr && opposite == nullptr) {
      throw std::invalid_argument("Portal has no links");
    }

    auto portal =
        std::make_shared<Portal>(gctx, std::move(along), std::move(opposite));
    portal->surface().assignGeometryId(
        surfacePointers.at(jPortal.at(kSurfaceIdKey))->geometryId());
    return portal;
  };

  PortalPointerLookup portalPointers;
  for (const auto& [portalId, record] : portalRecords) {
    const auto inserted =
        portalPointers.emplace(portalId, decodePortal(record.payload));
    if (!inserted) {
      throw std::invalid_argument("Portal pointer reconstruction failed");
    }
  }

  auto logger =
      Acts::getDefaultLogger("navigationPolicyLogger", Acts::Logging::INFO);
  for (const auto& [volumeId, record] : volumeRecords) {
    auto* volume = volumePointers.find(volumeId);
    if (volume == nullptr) {
      throw std::invalid_argument("Volume pointer reconstruction failed");
    }

    for (const std::size_t portalId : record.portalIds) {
      volume->addPortal(portalPointers.at(portalId));
    }
    for (const std::size_t surfaceId : record.surfaceIds) {
      volume->addSurface(surfacePointers.at(surfaceId));
    }

    volume->setNavigationPolicy(navigationPolicyFromJson(
        gctx, record.navigationPolicy, *volume, *logger));
  }

  auto root = std::move(volumeStorage.at(rootVolumeId));
  if (root == nullptr) {
    throw std::invalid_argument("Root volume reconstruction failed");
  }

  return std::shared_ptr<TrackingVolume>(std::move(root));
}

std::shared_ptr<Acts::TrackingGeometry>
Acts::TrackingGeometryJsonConverter::fromJson(const GeometryContext& gctx,
                                              const nlohmann::json& encoded,
                                              const Options& options) const {
  auto world = trackingVolumeFromJson(gctx, encoded, options);

  GeometryIdentifier::Value nextVolumeId = 1u;
  ensureIdentifiers(*world, nextVolumeId);

  return std::make_shared<TrackingGeometry>(
      world, nullptr, GeometryIdentifierHook{}, getDummyLogger(), false);
}
