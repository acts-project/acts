// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test_suite.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"

#include <filesystem>
#include <fstream>
#include <memory>

#include <Eigen/Geometry>
#include <nlohmann/json.hpp>

using namespace Acts;

const std::vector<std::shared_ptr<Surface>> surfaces = []() {
  std::vector<std::shared_ptr<Surface>> v;

  for (int i = 0; i < 3; ++i) {
    auto bounds = std::make_shared<RectangleBounds>(1.0, 1.0);
    Transform3 transform = Transform3::Identity();
    transform.translate(Vector3::Random());
    Vector3 randomAngles = Vector3::Random();
    SquareMatrix3 rotMatrix;
    rotMatrix = Eigen::AngleAxis(randomAngles[0], Vector3::UnitX()) *
                Eigen::AngleAxis(randomAngles[1], Vector3::UnitY()) *
                Eigen::AngleAxis(randomAngles[2], Vector3::UnitZ());
    transform.rotate(rotMatrix);
    v.push_back(Surface::makeShared<PlaneSurface>(transform, bounds));
  }

  return v;
}();

const std::string filename = "json_surfaces_reader_tests.json";

struct FileFixture {
  FileFixture() {
    nlohmann::json js = nlohmann::json::array();

    for (const auto &s : surfaces) {
      js.push_back(SurfaceJsonConverter::toJson({}, *s));
    }

    nlohmann::json j;
    j["foo"] = js;

    std::ofstream f(filename);
    f << j.dump(2);
  }

  ~FileFixture() { std::filesystem::remove(filename); }
};

FileFixture fileFixture;

namespace ActsTests {

const GeometryContext gctx{};

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(surface_reading_test) {
  auto readBackSurfaces = JsonSurfacesReader::readVector({filename, {"foo"}});

  BOOST_REQUIRE_EQUAL(surfaces.size(), readBackSurfaces.size());
  for (auto [refSurface, surface] : zip(surfaces, readBackSurfaces)) {
    BOOST_CHECK(refSurface->localToGlobalTransform(gctx).isApprox(
        surface->localToGlobalTransform(gctx), 1.e-4));
    BOOST_CHECK(
        refSurface->center(gctx).isApprox(surface->center(gctx), 1.e-4));
    BOOST_CHECK_EQUAL(refSurface->type(), surface->type());
    BOOST_CHECK_EQUAL(refSurface->bounds().type(), surface->bounds().type());
  }
}

BOOST_AUTO_TEST_CASE(json_detelement_reading_test) {
  auto readBackDetElements =
      JsonSurfacesReader::readDetectorElements({filename, {"foo"}}, 1.0);

  BOOST_REQUIRE_EQUAL(surfaces.size(), readBackDetElements.size());
  for (auto [refSurface, detElement] : zip(surfaces, readBackDetElements)) {
    auto surface = &detElement->surface();
    BOOST_CHECK(refSurface->localToGlobalTransform(gctx).isApprox(
        surface->localToGlobalTransform(gctx), 1.e-4));
    BOOST_CHECK(
        refSurface->center(gctx).isApprox(surface->center(gctx), 1.e-4));
    BOOST_CHECK_EQUAL(refSurface->type(), surface->type());
    BOOST_CHECK_EQUAL(refSurface->bounds().type(), surface->bounds().type());
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
