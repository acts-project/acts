// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/GeoModel/GeoModelModuleSplitter.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"

#include <GeoModelKernel/GeoBox.h>
#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>

using namespace Acts;

BOOST_AUTO_TEST_SUITE(GeoModelPlugin)

const GeometryContext gctx;

const auto box = GeoIntrusivePtr(new GeoBox(100, 200, 2));
const auto log = GeoIntrusivePtr(new GeoLogVol("LogVolume", box, nullptr));
const auto vol = GeoIntrusivePtr(new GeoFullPhysVol(log));

const std::array<double, AnnulusBounds::eSize> annulusParams{
    /* rmin      */ 384.0,
    /* rmax      */ 498,
    /* phiRelMin */ -0.07,
    /* phiRelMax */ 0.11,
    /* avgPhi    */ 0.0,
    /* shift_x   */ 5.0,
    /* shift_y   */ 8.0};

auto makeDetElement() {
  auto bounds = std::make_shared<AnnulusBounds>(annulusParams);

  return GeoModelDetectorElement::createDetectorElement<DiscSurface,
                                                        AnnulusBounds>(
      vol, bounds, Transform3::Identity(), 0.5);
}

BOOST_AUTO_TEST_CASE(ModuleSplitterTest_empty) {
  auto detEl = makeDetElement();
  const std::map<std::string, std::vector<double>> patterns{};

  GeoModelModuleSplitter splitter(patterns);

  BOOST_CHECK(splitter.split(detEl, gctx) == std::nullopt);
}

BOOST_AUTO_TEST_CASE(ModuleSplitterTest_non_annulus) {
  auto bounds = std::make_shared<RadialBounds>(1.0, 10.0);

  auto detEl =
      GeoModelDetectorElement::createDetectorElement<DiscSurface, RadialBounds>(
          vol, bounds, Transform3::Identity(), 0.5);

  const std::map<std::string, std::vector<double>> patterns{
      {"A", {1.0, 5.0, 10.0}}};
  GeoModelModuleSplitter splitter(patterns);

  BOOST_CHECK(splitter.split(detEl, gctx) == std::nullopt);
}

BOOST_AUTO_TEST_CASE(ModuleSplitterTest_split) {
  auto detEl = makeDetElement();
  auto midR = 0.5 * (annulusParams[1] - annulusParams[0]);
  const std::map<std::string, std::vector<double>> patterns{
      {"A", {annulusParams[0], midR, annulusParams[1]}}};

  double tolerance = 1.e-8;
  GeoModelModuleSplitter splitter(patterns, tolerance);

  auto res = splitter.split(detEl, gctx);
  BOOST_REQUIRE(res.has_value());
  BOOST_REQUIRE(res->size() == 2);

  BOOST_CHECK_CLOSE(res->at(0)->surface().bounds().values()[0],
                    annulusParams[0], tolerance);
  BOOST_CHECK_CLOSE(res->at(0)->surface().bounds().values()[1], midR,
                    tolerance);

  BOOST_CHECK_CLOSE(res->at(1)->surface().bounds().values()[0], midR,
                    tolerance);
  BOOST_CHECK_CLOSE(res->at(1)->surface().bounds().values()[1],
                    annulusParams[1], tolerance);
}

BOOST_AUTO_TEST_SUITE_END()
