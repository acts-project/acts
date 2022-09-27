// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/PortalSvgConverter.hpp"

#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"

#include <any>

namespace {

/// Helper method to make a proto link implementation
///
/// @param portalOptions are the options containin portal draw length and volume colour map
/// @param position is the start position of the link
/// @param direciton is the direction of the link
/// @param dVolume the volume this link is directing to
///
/// @return a protoLink
Acts::Svg::ProtoLink makeProtoLink(
    const Acts::Svg::PortalConverter::Options& portalOptions,
    const Acts::Vector3& position, const Acts::Vector3& direction,
    const Acts::Experimental::DetectorVolume* dVolume) {
  Acts::Svg::ProtoLink pLink;
  Acts::Vector3 end3 = position + portalOptions.linkLength * direction;
  pLink._start = {position.x(), position.y(), position.z()};
  pLink._end = {end3.x(), end3.y(), end3.z()};
  auto linkIndexCandidate = portalOptions.volumeIndices.find(dVolume);
  if (linkIndexCandidate != portalOptions.volumeIndices.end()) {
    pLink._link_index = linkIndexCandidate->second;
  }
  return pLink;
};

/// Helper method to convert a multi link
///
/// @param portalOptions are the options containin portal draw length and volume colour map
/// @param multiLink is the mulit link to be converted
/// @param surface is the portal surface
/// @param refPosition is the reference position which might be needed to decode
/// @param backward boolean to tell if its a backward link
/// @param gctx is the geometry context
///
/// @return it will return the proto links
std::vector<Acts::Svg::ProtoLink> convertMultiLink(
    const Acts::Svg::PortalConverter::Options& portalOptions,
    const Acts::Experimental::detail::MultiLink1DImpl& multiLink,
    const Acts::Surface& surface, const Acts::Vector3& refPosition,
    bool backward, const Acts::GeometryContext& gctx) noexcept(false) {
  // The return links
  std::vector<Acts::Svg::ProtoLink> pLinks;
  int sign = backward ? -1 : 1;
  // Generate the proto-links of the multi-link
  for (auto [il, v] : Acts::enumerate(multiLink.dVolumes)) {
    Acts::Vector3 position = refPosition;
    Acts::ActsScalar refC =
        0.5 * (multiLink.cBoundaries[il + 1u] + multiLink.cBoundaries[il]);
    if (multiLink.bValue == Acts::binR) {
      Acts::ActsScalar phi = Acts::VectorHelpers::phi(refPosition);
      position = Acts::Vector3(refC * std::cos(phi), refC * std::sin(phi),
                               refPosition.z());
    } else if (multiLink.bValue == Acts::binZ) {
      position[2] = refC;
    } else if (multiLink.bValue == Acts::binPhi) {
      Acts::ActsScalar r = Acts::VectorHelpers::perp(refPosition);
      position = Acts::Vector3(r * std::cos(refC), r * std::sin(refC),
                               refPosition.z());
    } else {
      throw std::invalid_argument("convertMultiLink: incorrect binning.");
    }
    Acts::Vector3 direction = surface.normal(gctx, position);
    pLinks.push_back(makeProtoLink(portalOptions, position,
                                   Acts::Vector3(sign * direction), v));
  }
  return pLinks;
}

}  // namespace

Acts::Svg::ProtoPortal Acts::Svg::PortalConverter::convert(
    const GeometryContext& gctx, const Experimental::Portal& portal,
    const PortalConverter::Options& portalOptions) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("PortalSvgConverter", portalOptions.logLevel));

  ProtoPortal pPortal;
  // First convert the surface
  pPortal._surface = SurfaceConverter::convert(gctx, portal.surface(),
                                               portalOptions.surfaceOptions);

  // Reference point and direction
  Vector3 rPos(0., 0., 0);
  Vector3 rDir(0., 0., 1);
  const auto& surface = portal.surface();
  const auto surfaceTransform = portal.surface().transform(gctx);
  const auto surfaceRotation = surfaceTransform.rotation().eval();
  const auto surfaceTranslation = surfaceTransform.translation().eval();
  const auto surfaceType = surface.bounds().type();
  const auto& boundValues = surface.bounds().values();
  switch (surfaceType) {
    case SurfaceBounds::eCylinder: {
      // Get phi
      ActsScalar r = boundValues[0u];
      ActsScalar aphi = boundValues[3u];
      rPos = Vector3(r * std::cos(aphi), r * std::sin(aphi),
                     surfaceTranslation.z());
    } break;
    case SurfaceBounds::eDisc: {
      // Get phi
      ActsScalar r = 0.5 * (boundValues[0u] + boundValues[1u]);
      ActsScalar aphi = boundValues[3u];
      rPos = Vector3(r * std::cos(aphi), r * std::sin(aphi),
                     surfaceTranslation.z());
    } break;
    default: {
      rPos = surfaceTranslation;
    } break;
  }
  rDir = surface.normal(gctx, rPos);

  // Now convert the link objects
  const auto [backwardLink, forwardLink] = portal.volumeLinks();

  /// @brief convert the managed link to a proto link for displayeing
  ///
  /// @param mvl the managed detector volume link (link, implemenation)
  /// @param backward is a flag indicating whether this is a backward link
  ///
  /// @return a vector of proto links
  auto convertLink = [&](const Experimental::ManagedDetectorVolumeLink& mvl,
                         bool backward = false) -> std::vector<ProtoLink> {
    if (mvl.implementation != nullptr) {
      // Define the sign
      int sign = backward ? -1 : 1;
      // dynamic cast cascade now ...
      auto singleLink =
          dynamic_cast<Acts::Experimental::detail::SingleLinkImpl*>(
              mvl.implementation.get());
      if (singleLink != nullptr) {
        return {makeProtoLink(portalOptions, rPos, Vector3(sign * rDir),
                              singleLink->dVolume)};
      }
      // Check for a multi link here
      auto multiLink =
          dynamic_cast<Acts::Experimental::detail::MultiLink1DImpl*>(
              mvl.implementation.get());
      if (multiLink != nullptr) {
        return convertMultiLink(portalOptions, *multiLink, surface, rPos,
                                backward, gctx);
      }
      // Check for a composed link
    }
    return {};
  };

  auto bwProtoLinks = convertLink(backwardLink, true);
  auto fwProtoLinks = convertLink(forwardLink);

  pPortal._volume_links.insert(pPortal._volume_links.end(),
                               bwProtoLinks.begin(), bwProtoLinks.end());
  pPortal._volume_links.insert(pPortal._volume_links.end(),
                               fwProtoLinks.begin(), fwProtoLinks.end());

  // Return the proto Portal
  return pPortal;
}
