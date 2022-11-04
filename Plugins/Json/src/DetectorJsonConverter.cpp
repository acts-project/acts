// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DetectorJsonConverter.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Experimental/Detector.hpp"
#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Experimental/detail/DetectorVolumeFinders.hpp"
#include "Acts/Experimental/detail/NavigationStateUpdators.hpp"
#include "Acts/Experimental/detail/SurfaceGridGenerator.hpp"
#include "Acts/Experimental/detail/SurfaceGridHelper.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/DetrayJsonHelper.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/PortalJsonConverter.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"
#include "Acts/Plugins/Json/VolumeBoundsJsonConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <iostream>
#include <memory>
#include <vector>

namespace {

/// @brief  A portal assinger class that acting as a portal generator in the
/// Detector geometry building
///
/// @note if no indices are given, all portals are assigned
struct PortalAssigner {
  const std::vector<std::shared_ptr<Acts::Experimental::Portal>>& portals;
  const std::vector<std::size_t>& indices;

  /// @brief  Constructor that takes references in order to avoid copying
  PortalAssigner(
      const std::vector<std::shared_ptr<Acts::Experimental::Portal>>& ps,
      const std::vector<std::size_t>& is = {})
      : portals(ps), indices(is) {}

  // The portal provider function just assigns
  std::vector<std::shared_ptr<Acts::Experimental::Portal>> assign(
      [[maybe_unused]] const Acts::Transform3& dTransform,
      [[maybe_unused]] const Acts::VolumeBounds& dBounds,
      [[maybe_unused]] std::shared_ptr<Acts::Experimental::DetectorVolume>
          dVolume) const {
    std::vector<std::shared_ptr<Acts::Experimental::Portal>> rportals;
    if (not indices.empty()) {
      rportals.reserve(indices.size());
      for (auto idx : indices) {
        rportals.push_back(portals[idx]);
      }
      return rportals;
    }
    return portals;
  }
};

/// @brief  Create a detector volume from json
///
/// @param j the json object
/// @param portals [in,out] the already created portals (if not detray writer)
/// @param jPortals [in,out]the updated json description of portals (for detray mode)
/// @param logLevel the log level
///
/// @note in detray mode, the portals are filled with the newly created ones
///
/// @return
std::shared_ptr<Acts::Experimental::DetectorVolume> detectorVolumeFromJson(
    const nlohmann::json& j,
    std::vector<std::shared_ptr<Acts::Experimental::Portal>>& portals,
    std::vector<nlohmann::json>& jPortals, Acts::Logging::Level logLevel) {
  ACTS_LOCAL_LOGGER(Acts::getDefaultLogger("JSON I/O:", logLevel));

  // Get the name of the volume
  std::string name = j["name"];
  ACTS_DEBUG("Reading volume " << name);

  // Create the bounds and the transform
  Acts::Transform3 transform;
  Acts::from_json(j["transform"], transform);
  auto volumeBounds = Acts::unqiueVolumeBoundsFromJson(j["bounds"]);

  // Create all the surfaces
  std::vector<std::shared_ptr<Acts::Surface>> surfaces;
  if (j.find("surfaces") != j.end()) {
    auto jSurfaces = j["surfaces"];
    for (const auto& jSurface : jSurfaces) {
      surfaces.push_back(Acts::surfaceFromJson(jSurface));
    }
  }

  std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes = {};

  // The volume portals and volume portal indices
  std::vector<std::shared_ptr<Acts::Experimental::Portal>> vPortals;
  std::vector<std::size_t> vPortalIndices = {};

  // Get the portal indices - in ACTS mode
  if (j.find("portalIndices") != j.end()) {
    ACTS_DEBUG("Portal indices present - reading in ACTS mode.");
    vPortals = portals;
    std::vector<std::size_t> aPortalIndices = j["portalIndices"];
    vPortalIndices = aPortalIndices;
  } else if (j.find("portals") != j.end()) {
    ACTS_DEBUG("Portals present - reading in detray mode.");
    const auto& vjPortals = j["portals"];
    for (const auto& jPortal : vjPortals) {
      if (not jPortal.empty()) {
        vPortals.push_back(Acts::Experimental::portalFromJson(jPortal));
      }
    }
    ACTS_VERBOSE("- Read " << vPortals.size()
                           << " volume portals (detray mode).")
    portals.insert(portals.end(), vPortals.begin(), vPortals.end());
    jPortals.insert(jPortals.end(), vjPortals.begin(), vjPortals.end());
  } else {
    ACTS_ERROR("Could not determine mode (ACTS/detray), bailing out.");
    return nullptr;
  }

  PortalAssigner pa(vPortals, vPortalIndices);
  Acts::Delegate<std::vector<std::shared_ptr<Acts::Experimental::Portal>>(
      const Acts::Transform3&, const Acts::VolumeBounds&,
      std::shared_ptr<Acts::Experimental::DetectorVolume>)>
      pAssigner;
  pAssigner.connect<&PortalAssigner::assign>(&pa);

  // Sort out local navigation
  std::shared_ptr<Acts::Experimental::DetectorVolume> volume =
      Acts::Experimental::DetectorVolumeFactory::construct(
          pAssigner, Acts::GeometryContext(), name, transform,
          std::move(volumeBounds), Acts::Experimental::detail::allPortals());

  if (j.find("surfaceNavigation") != j.end()) {
    // The surface navigation
    nlohmann::json jsNavigation = j["surfaceNavigation"];
    // @ hack for now just to get it going
    if (jsNavigation["type"] == "grid") {
      std::vector<std::size_t> globalBins = jsNavigation["globalBins"];
      std::vector<std::vector<std::size_t>> entries = jsNavigation["entries"];
      std::array<Acts::BinningValue, 2u> bValues;
      std::vector<Acts::ActsScalar> mins;
      std::vector<Acts::ActsScalar> maxs;
      std::vector<std::size_t> bins;

      auto jAxes = jsNavigation["axes"];
      for (auto [ij, jAxis] : Acts::enumerate(jAxes)) {
        auto ja = jAxis["axis"];
        std::array<Acts::ActsScalar, 2u> range = ja["range"];
        mins.push_back(range[0u]);
        maxs.push_back(range[1u]);
        bins.push_back(ja["bins"]);
        std::string bValStr = jAxis["binningValue"];
        if (bValStr == "binPhi") {
          bValues[ij] = Acts::binPhi;
        } else if (bValStr == "binR") {
          bValues[ij] = Acts::binR;
        } else {
          bValues[ij] = Acts::binZ;
        }
      }
      Acts::Experimental::detail::AxisEqB rzA(mins[0u], maxs[0u], bins[0u]);
      Acts::Experimental::detail::AxisEqC phiA(mins[1u], maxs[1u], bins[1u]);
      Acts::Experimental::detail::GridEqBEqC rzPhiGrid(
          std::make_tuple(std::move(rzA), std::move(phiA)));

      for (auto [i, igb] : Acts::enumerate(globalBins)) {
        rzPhiGrid.at(igb) = entries[i];
      }

      // Create the surface navigation and attach
      Acts::Experimental::ManagedNavigationStateUpdator managedUpdator;

      // Create grid surface and all portals attacher
      auto gSurfaceAttacher = Acts::Experimental::detail::GridSurfaceAttacher<
          Acts::Experimental::detail::GridEqBEqC>(
          {std::move(rzPhiGrid), bValues});

      auto portalAttacher = Acts::Experimental::detail::AllPortalsAttacher{};

      // Combine to a navigation state updator
      using NavigationUpdatorImpl =
          Acts::Experimental::detail::NavigationStateUpdatorImpl<
              decltype(gSurfaceAttacher), decltype(portalAttacher)>;

      auto navUpdator = std::make_shared<NavigationUpdatorImpl>(
          std::make_tuple(std::move(gSurfaceAttacher), portalAttacher));

      Acts::Experimental::NavigationStateUpdator nStateUpdator;
      nStateUpdator.connect<&NavigationUpdatorImpl::update>(navUpdator.get());
      managedUpdator.delegate = std::move(nStateUpdator);
      managedUpdator.implementation = navUpdator;

      // Update the navigaiton state updator
      volume->updateNavigationStateUpator(std::move(managedUpdator), surfaces);
    }
  }
  return volume;
}

};  // namespace

nlohmann::json Acts::Experimental::toJson(
    const Acts::Experimental::Detector& detector,
    const Acts::GeometryContext& gctx, bool detray,
    Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("JSON I/O:", logLevel));

  // The return object
  nlohmann::json jDetector;
  jDetector["name"] = detector.name();
  // Conatined volumes
  auto volumes = detector.volumes();

  // Make a portal map first
  std::map<const Acts::Experimental::Portal*, unsigned int> portalIndexMap;
  if (not detray) {
    unsigned int pIndex = 0u;
    for (const auto& v : volumes) {
      for (const auto& p : v->portals()) {
        if (portalIndexMap.find(p) == portalIndexMap.end()) {
          portalIndexMap[p] = pIndex++;
        }
      }
    }
    nlohmann::json jPortals;
    for (auto [p, pindex] : portalIndexMap) {
      jPortals.push_back(toJson(*p, volumes, gctx));
    }
    jDetector["portals"] = jPortals;
  }
  // Get the list of volumes
  nlohmann::json jVolumes;
  for (auto [iv, v] : enumerate(volumes)) {
    nlohmann::json jVolume;
    jVolume["name"] = v->name();
    jVolume["bounds"] = v->volumeBounds();
    nlohmann::json jTransform;
    Acts::to_json(jTransform, v->transform(gctx));
    jVolume["transform"] = jTransform;
    // Deal with portals --------------------------
    nlohmann::json jPortals;
    // for detray make unique portals
    if (detray) {
      // We need to clip and adjust the portals
      auto clippedPortals = DetrayJsonHelper::clipAndPickPortals(gctx, *v);
      // Write the portals per volume
      for (auto [ip, p] : enumerate(clippedPortals)) {
        jPortals.push_back(toJson(*p, volumes, gctx));
      }
      jVolume["portals"] = jPortals;
    } else {
      // Write the portal indices are they might be shared
      for (const auto& p : v->portals()) {
        if (portalIndexMap.find(p) != portalIndexMap.end()) {
          jPortals.push_back(portalIndexMap.find(p)->second);
        }
      }
      jVolume["portalIndices"] = jPortals;
    }
    // Deal with surfaces --------------------------
    if (not v->surfaces().empty()) {
      nlohmann::json jSurfaces;
      for (auto [is, s] : enumerate(v->surfaces())) {
        nlohmann::json jSurface;
        toJson(jSurface, *s, gctx);
        jSurfaces.push_back(jSurface);
      }
      jVolume["surfaces"] = jSurfaces;
    }

    // Write the managed navigation state updator
    // Get access to the navigation state updator
    auto nStateUpdator = v->navigationStateUpdator();
    auto nStateUpdatorImpl = nStateUpdator.implementation;
    if (nStateUpdatorImpl != nullptr) {
      // Get the grid parameters
      auto gParameters = Experimental::detail::extractGridParameters(
          *nStateUpdatorImpl, Experimental::detail::s_gridUpdatorTempaltes);
      if (gParameters.axes.size() == 2u) {
        nlohmann::json jsNavigation;
        jsNavigation["type"] = "grid";
        jsNavigation["entries"] = gParameters.entries;
        jsNavigation["globalBins"] = gParameters.globalBins;
        nlohmann::json jsAxes;
        for (const auto& ga : gParameters.axes) {
          nlohmann::json ajs;
          ajs["binningValue"] = binningValueNames()[ga.bValue];
          ajs["axis"] = *ga.raw;
          jsAxes.push_back(ajs);
        }
        jsNavigation["axes"] = jsAxes;
        jVolume["surfaceNavigation"] = jsNavigation;
      }
    }
    // Add to the list of volumes
    jVolumes.push_back(jVolume);
  }
  // Add to the detector
  jDetector["volumes"] = jVolumes;
  // Return the detector
  return jDetector;
}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::detectorFromJson(const nlohmann::json& j,
                                     Acts::Logging::Level logLevel) {
  // The local logger
  ACTS_LOCAL_LOGGER(getDefaultLogger("JSON I/O:", logLevel));

  if (j.find("detector") != j.end()) {
    auto jDetector = j["detector"];
    ACTS_INFO("Reading detector with name " << jDetector["name"]);

    std::vector<std::shared_ptr<Acts::Experimental::Portal>> portals;
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes;

    std::vector<nlohmann::json> jPortals;
    // Portal reading section
    if (jDetector.find("portals") != jDetector.end()) {
      ACTS_INFO("Portal section present - reading in ACTS mode.");
      std::vector<nlohmann::json> njPortals = jDetector["portals"];
      jPortals = njPortals;  // nlohmann needs explicit type for reading
      for (const auto& jPortal : jPortals) {
        if (not jPortal.empty()) {
          portals.push_back(portalFromJson(jPortal));
        }
      }
    } else {
      ACTS_INFO("Portal section NOT present - reading in detray mode.");
    }

    // Volume reading section
    if (jDetector.find("volumes") != jDetector.end()) {
      auto jVolumes = jDetector["volumes"];
      for (const auto& jVolume : jVolumes) {
        volumes.push_back(
            detectorVolumeFromJson(jVolume, portals, jPortals, logLevel));
      }
      ACTS_INFO("- Read " << volumes.size() << " volumes.");
      ACTS_INFO("- Read " << portals.size() << " portals.");
    }

    // Portals and volumes exist - all the links can be set
    if (not portals.empty() and not jPortals.empty() and not volumes.empty()) {
      for (auto [ip, jPortal] : enumerate(jPortals)) {
        auto& portal = portals[ip];
        attachVolumeLinkFromJson(jPortal, portal, volumes, logLevel);
      }
    }

    // A detector construction that should world
    return Detector::makeShared(jDetector["name"], volumes,
                                detail::tryAllVolumes());
  }

  return nullptr;
}
