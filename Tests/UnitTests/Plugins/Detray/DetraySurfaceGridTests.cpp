// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Detray/DetrayJsonConverter.hpp"
#include "ActsPlugins/Detray/DetrayPayloadConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceGridPolicy.hpp"
#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"
#include "ActsTests/CommonHelpers/CylindricalTrackingGeometry.hpp"

#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <detray/geometry/identifier.hpp>
#include <detray/io/frontend/payloads.hpp>

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace ActsPlugins;

namespace ActsTests {

namespace {

auto gctx = GeometryContext::dangerouslyDefaultConstruct();

std::shared_ptr<TrackingGeometry> buildGeometry(
    std::optional<DetraySurfaceGridConfig> policyConfig = std::nullopt) {
  // The builder owns the detector elements, keep it alive for the test
  static std::vector<std::unique_ptr<CylindricalTrackingGeometry>> builders;
  auto& cGeometry = *builders.emplace_back(
      std::make_unique<CylindricalTrackingGeometry>(gctx, true));
  if (policyConfig.has_value()) {
    cGeometry.extendLayerNavigationPolicy =
        [cfg = *policyConfig](NavigationPolicyFactory&& factory) {
          return std::move(factory).add<DetraySurfaceGridPolicy>(cfg);
        };
  }
  return cGeometry(*getDefaultLogger("Geo", Logging::WARNING));
}

/// Layer name -> bin index -> surface indices
using GridContents =
    std::map<std::string,
             std::map<std::pair<unsigned, unsigned>, std::set<std::size_t>>>;

struct Converted {
  GridContents grids;
  /// Layer name -> sensitive surface indices
  std::map<std::string, std::set<std::size_t>> sensitives;
};

Converted convert(const TrackingGeometry& geometry,
                  std::optional<DetraySurfaceGridConfig> converterFill) {
  DetrayPayloadConverter::Config cfg;
  cfg.surfaceGridFill = converterFill;
  geometry.apply([&cfg](const TrackingVolume& volume) {
    if (volume.volumeName() == "Beampipe") {
      cfg.beampipeVolume = &volume;
    }
  });
  DetrayPayloadConverter converter(cfg,
                                   getDefaultLogger("Cnv", Logging::WARNING));
  auto payloads = converter.convertTrackingGeometry(gctx, geometry);

  std::map<std::size_t, std::string> names;
  Converted out;
  for (const auto& volume : payloads.geometry.volumes) {
    names[volume.index.link] = volume.name;
    for (const auto& [i, surface] : enumerate(volume.surfaces)) {
      if (surface.type == detray::surface_id::e_sensitive) {
        out.sensitives[volume.name].insert(i);
      }
    }
  }
  for (const auto& [volumeIndex, grids] : payloads.surface_grids->grids) {
    BOOST_REQUIRE_EQUAL(grids.size(), 1u);
    auto& contents = out.grids[names.at(volumeIndex)];
    for (const auto& bin : grids.front().bins) {
      contents[{bin.loc_index.at(0), bin.loc_index.at(1)}].insert(
          bin.content.begin(), bin.content.end());
    }
  }
  return out;
}

std::size_t entries(const GridContents& grids) {
  std::size_t n = 0;
  for (const auto& [name, bins] : grids) {
    for (const auto& [bin, content] : bins) {
      n += content.size();
    }
  }
  return n;
}

DetraySurfaceGridConfig straightFromOrigin() {
  DetraySurfaceGridConfig cfg;
  cfg.vertexZHalfRange = 0.;
  cfg.nVertexZ = 1;
  cfg.minPt = 0.;
  return cfg;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(DetraySuite)

BOOST_AUTO_TEST_CASE(DetraySurfaceGridFillReachesEverySensitive) {
  auto geometry = buildGeometry();
  const Converted filled = convert(*geometry, DetraySurfaceGridConfig{});
  const Converted copied = convert(*geometry, std::nullopt);

  BOOST_CHECK_EQUAL(filled.grids.size(), 4u);
  for (const auto& [name, bins] : filled.grids) {
    std::set<std::size_t> seen;
    for (const auto& [bin, content] : bins) {
      seen.insert(content.begin(), content.end());
    }
    BOOST_CHECK(seen == filled.sensitives.at(name));
    // Same binning as the copied ACTS grid, the explicit bins are kept
    BOOST_CHECK_EQUAL(bins.size(), copied.grids.at(name).size());
  }

  // Filled for detray's lookup, not a copy of the layer footprint
  BOOST_CHECK(filled.grids != copied.grids);
}

BOOST_AUTO_TEST_CASE(DetraySurfaceGridFillFamilyWidensBins) {
  auto geometry = buildGeometry();
  const Converted narrow = convert(*geometry, straightFromOrigin());
  const Converted wide = convert(*geometry, DetraySurfaceGridConfig{});

  // A wider track family only adds surfaces to bins
  for (const auto& [name, bins] : narrow.grids) {
    for (const auto& [bin, content] : bins) {
      const auto& wider = wide.grids.at(name).at(bin);
      BOOST_CHECK(std::ranges::includes(wider, content));
    }
  }
  BOOST_CHECK_GT(entries(wide.grids), entries(narrow.grids));
}

BOOST_AUTO_TEST_CASE(DetraySurfaceGridPolicyOverridesConverter) {
  // The policy in the geometry wins over the converter default
  auto withPolicy = buildGeometry(straightFromOrigin());
  const Converted fromPolicy = convert(*withPolicy, DetraySurfaceGridConfig{});

  auto plain = buildGeometry();
  const Converted fromConverter = convert(*plain, straightFromOrigin());

  BOOST_CHECK(fromPolicy.grids == fromConverter.grids);

  // ... also when the converter would otherwise copy the ACTS bins
  const Converted fromPolicyOnly = convert(*withPolicy, std::nullopt);
  BOOST_CHECK(fromPolicyOnly.grids == fromConverter.grids);
}

BOOST_AUTO_TEST_CASE(DetraySurfaceGridPolicyJsonRoundTrip) {
  DetraySurfaceGridConfig cfg = straightFromOrigin();
  cfg.bField = 3_T;
  cfg.numberOfBinsFactor = 1.5;
  auto geometry = buildGeometry(cfg);

  // Without registration the writer does not know the policy
  TrackingGeometryJsonConverter plainWriter(
      TrackingGeometryJsonConverter::Config::defaultConfig());
  BOOST_CHECK_THROW(plainWriter.toJson(gctx, *geometry), std::exception);

  auto jsonCfg = TrackingGeometryJsonConverter::Config::defaultConfig();
  addDetrayNavigationPolicies(jsonCfg);
  TrackingGeometryJsonConverter converter(jsonCfg);
  auto readBack = converter.fromJson(gctx, converter.toJson(gctx, *geometry));

  std::size_t nPolicies = 0;
  readBack->apply([&](const TrackingVolume& volume) {
    const INavigationPolicy* policy = volume.navigationPolicy();
    if (policy == nullptr) {
      return;
    }
    policy->visit([&](const INavigationPolicy& p) {
      if (const auto* grid = dynamic_cast<const DetraySurfaceGridPolicy*>(&p);
          grid != nullptr) {
        ++nPolicies;
        BOOST_CHECK_EQUAL(grid->config().vertexZHalfRange,
                          cfg.vertexZHalfRange);
        BOOST_CHECK_EQUAL(grid->config().minPt, cfg.minPt);
        BOOST_CHECK_EQUAL(grid->config().bField, cfg.bField);
        BOOST_CHECK_EQUAL(grid->config().nVertexZ, cfg.nVertexZ);
        BOOST_CHECK_EQUAL(grid->config().nCurvature, cfg.nCurvature);
        BOOST_CHECK(grid->config().numberOfBinsFactor ==
                    cfg.numberOfBinsFactor);
      }
    });
  });
  BOOST_CHECK_EQUAL(nPolicies, 4u);

  // The read-back geometry converts to the same grids
  BOOST_CHECK(convert(*readBack, std::nullopt).grids ==
              convert(*geometry, std::nullopt).grids);
}

BOOST_AUTO_TEST_CASE(DetraySurfaceGridLookupRadius) {
  auto geometry = buildGeometry();

  auto radii = [&](bool onLayer) {
    DetrayPayloadConverter::Config cfg;
    DetraySurfaceGridConfig fill;
    fill.lookupOnLayer = onLayer;
    cfg.surfaceGridFill = fill;
    geometry->apply([&cfg](const TrackingVolume& volume) {
      if (volume.volumeName() == "Beampipe") {
        cfg.beampipeVolume = &volume;
      }
    });
    DetrayPayloadConverter converter(cfg,
                                     getDefaultLogger("Cnv", Logging::WARNING));
    auto payloads = converter.convertTrackingGeometry(gctx, *geometry);
    std::vector<std::optional<double>> out;
    for (const auto& [volume, grids] : payloads.surface_grids->grids) {
      out.push_back(grids.front().reference_radius);
    }
    return out;
  };

  // On the layer: the mean module radius, inside the layer volume
  const auto onLayer = radii(true);
  BOOST_CHECK_EQUAL(onLayer.size(), 4u);
  for (const auto& r : onLayer) {
    BOOST_REQUIRE(r.has_value());
    BOOST_CHECK_GT(*r, 0.);
  }
  // Otherwise detray derives it from the portals
  for (const auto& r : radii(false)) {
    BOOST_CHECK(!r.has_value());
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
