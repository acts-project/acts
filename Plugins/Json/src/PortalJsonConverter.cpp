// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/PortalJsonConverter.hpp"

#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"

nlohmann::json Acts::Experimental::toJson(
    const Portal& portal, const std::vector<const DetectorVolume*>& volumes,
    const GeometryContext& gctx, bool detray, Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("JSON I/O:", logLevel));

  ACTS_VERBOSE("- Writing a portal in " << (detray ? "detray" : "ACTS")
                                        << " mode.");

  // The overall return object
  nlohmann::json jPortal;

  nlohmann::json jSurface;
  toJson(jSurface, portal.surface(), gctx);
  jPortal["surface"] = jSurface;

  // The Portal link
  nlohmann::json jLinks;
  auto volumeLinks = portal.volumeLinks();
  for (auto [il, vl] : enumerate(volumeLinks)) {
    if (vl.implementation != nullptr) {
      jLinks.push_back(toJson(*vl.implementation, volumes));
    } else {
      nlohmann::json eow;
      eow["type"] = "endOfWorld";
      jLinks.push_back(eow);
    }
  }
  jPortal["volumeLinks"] = jLinks;
  // Return the full json object
  return jPortal;
}

nlohmann::json Acts::Experimental::toJson(
    const IDelegateImpl& delegate,
    const std::vector<const DetectorVolume*>& volumes,
    Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("JSON I/O:", logLevel));

  auto findVolume = [&](const DetectorVolume* volume) -> int {
    auto candidate = std::find(volumes.begin(), volumes.end(), volume);
    if (candidate != volumes.end()) {
      return std::distance(volumes.begin(), candidate);
    }
    return -1;
  };

  nlohmann::json jLink;
  // dynamic cast cascade
  const auto* delPtr = &delegate;
  // Check if single link present
  auto singleLink = dynamic_cast<const detail::SingleLinkImpl*>(delPtr);
  if (singleLink != nullptr) {
    int vIndex = findVolume(singleLink->dVolume);
    jLink["type"] = "single";
    jLink["link"] = vIndex;
    return jLink;
  }
  /// Check if multiple link  or transformed multilink is present
  ///
  /// @param mLink the multi link w/o transform
  /// @param transform an optional transform
  auto writeMultiLink1D = [&](const detail::MultiLink1DImpl& mLink,
                              const Transform3* transform = nullptr) -> void {
    jLink["type"] = (transform == nullptr) ? "multi1D" : "multi1D_transformed";
    jLink["binning"] = binningValueNames()[mLink.bValue];
    jLink["boundaries"] = mLink.cBoundaries;
    // Now find the volume links
    std::vector<unsigned int> vIndices = {};
    for (const auto& v : mLink.dVolumes) {
      vIndices.push_back(findVolume(v));
    }
    jLink["links"] = vIndices;
    if (transform != nullptr) {
      nlohmann::json jTransform;
      Acts::to_json(jTransform, *transform);
      jLink["transform"] = jTransform;
    }
  };

  // It is a multi link
  auto multiLink1D = dynamic_cast<const detail::MultiLink1DImpl*>(delPtr);
  if (multiLink1D != nullptr) {
    writeMultiLink1D(*multiLink1D);
    return jLink;
  }

  // It is a multi link with transfrom
  auto tMultiLink1D =
      dynamic_cast<const detail::TransformedMultiLink1DImpl*>(delPtr);
  if (tMultiLink1D != nullptr) {
    writeMultiLink1D(tMultiLink1D->multiLink, &(tMultiLink1D->transform));
  }

  return jLink;
}

std::shared_ptr<Acts::Experimental::Portal> Acts::Experimental::portalFromJson(
    const nlohmann::json& j, Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("JSON I/O:", logLevel));

  // Create the portal + surface
  auto jSurface = j["surface"];
  auto portal =
      Acts::Experimental::Portal::makeShared(surfaceFromJson(jSurface));

  return portal;
}

void Acts::Experimental::attachVolumeLinkFromJson(
    const nlohmann::json& j, std::shared_ptr<Portal>& portal,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("JSON I/O:", logLevel));

  auto jLinks = j["volumeLinks"];

  ACTS_VERBOSE("- Attaching a portal link(s) in "
               << (jLinks.size() == 1u ? "detray" : "ACTS") << " mode.");

  // Create the links
  for (auto [ij, jLink] : enumerate(jLinks)) {
    // The navigation direction
    Acts::NavigationDirection nDir = ij == 0
                                         ? Acts::NavigationDirection::Backward
                                         : Acts::NavigationDirection::Forward;

    std::string lType = jLink["type"];
    if (lType == "endOfWorld") {
      continue;
    } else if (lType == "single") {
      std::size_t link = jLink["link"];
      // Create a single link
      auto singleLinkImpl =
          std::make_shared<detail::SingleLinkImpl>(*volumes[link].get());
      DetectorVolumeLink singleLink;
      singleLink.connect<&detail::SingleLinkImpl::targetVolume>(
          singleLinkImpl.get());
      ManagedDetectorVolumeLink managedLink{std::move(singleLink),
                                            std::move(singleLinkImpl)};
      portal->updateVolumeLink(nDir, std::move(managedLink), {volumes[link]});
    } else if (lType == "multi1D") {
      // Let's create a mulit link
      std::vector<ActsScalar> boundaries = jLink["boundaries"];
      std::vector<std::size_t> links = jLink["links"];
      std::vector<std::shared_ptr<DetectorVolume>> attachedVolumes;
      for (const auto& l : links) {
        attachedVolumes.push_back(volumes[l]);
      }
      auto constVolumes = unpack_shared_const_vector(attachedVolumes);
      /// @brief  @TODO write a helper for that
      std::string binningValue = jLink["binning"];
      const auto bValueNames = binningValueNames();
      auto bvItr =
          std::find(bValueNames.begin(), bValueNames.end(), binningValue);
      if (bvItr != bValueNames.end()) {
        BinningValue bValue =
            s_binningValues[std::distance(bValueNames.begin(), bvItr)];

        auto multiLinkImpl = std::make_shared<detail::MultiLink1DImpl>(
            constVolumes, boundaries, bValue);
        // Create the delegate and its managed object
        DetectorVolumeLink multiLink;
        multiLink.connect<&detail::MultiLink1DImpl::targetVolume>(
            multiLinkImpl.get());
        ManagedDetectorVolumeLink managedLink =
            Acts::Experimental::ManagedDetectorVolumeLink{std::move(multiLink),
                                                          multiLinkImpl};
        portal->updateVolumeLink(nDir, std::move(managedLink), attachedVolumes);
      }
    } else if (lType == "multi1D_transformed") {
      std::cout << " not implemented yet" << std::endl;
    }
  }
}
