// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"
#include "ActsPlugins/DD4hep/DD4hepMaterialHelpers.hpp"

#include <numbers>

using namespace Acts;
using namespace ActsPlugins;

BOOST_AUTO_TEST_CASE(DD4hepProtoMaterialDefersRanges) {
  dd4hep::rec::VariantParameters params;
  // Construct variants in place to avoid Boost variant assignment warnings
  // with optimized GCC builds.
  params.variantParameters.emplace("material_binPhi", 8);
  params.variantParameters.emplace("material_binZ", 3);
  params.variantParameters.emplace("material_binR", 4);

  const auto cylinderMaterial = createProtoMaterial(
      params, "material", {{"binPhi", closed}, {"binZ", open}});
  const auto discMaterial = createProtoMaterial(
      params, "material", {{"binPhi", closed}, {"binR", open}});
  for (const auto& material : {cylinderMaterial, discMaterial}) {
    BOOST_CHECK(material->binning().isDeferred());
    for (const auto& axis : material->binning().axisSpecs()) {
      BOOST_CHECK(!axis.asEquidistant().min);
      BOOST_CHECK(!axis.asEquidistant().max);
    }
  }
  BOOST_CHECK(cylinderMaterial->binning().axisSpec(0).direction() ==
              AxisDirection::AxisRPhi);
  BOOST_CHECK(discMaterial->binning().axisSpec(0).direction() ==
              AxisDirection::AxisPhi);

  const auto cylinder =
      Surface::makeShared<CylinderSurface>(Transform3::Identity(), 30., 100.);
  const auto cylinderAxes =
      resolveMultiAxis(cylinderMaterial->binning(), *cylinder);
  BOOST_CHECK_EQUAL(cylinderAxes->getAxis(0).getNBins(), 8);
  BOOST_CHECK_CLOSE(cylinderAxes->getAxis(0).getMax(), 30. * std::numbers::pi,
                    1e-10);
  BOOST_CHECK_EQUAL(cylinderAxes->getAxis(1).getNBins(), 3);
  BOOST_CHECK_EQUAL(cylinderAxes->getAxis(1).getMin(), -100.);
  BOOST_CHECK_EQUAL(cylinderAxes->getAxis(1).getMax(), 100.);

  const auto disc = Surface::makeShared<DiscSurface>(
      Transform3::Identity(), std::make_shared<RadialBounds>(10., 50.));
  const auto discAxes = resolveMultiAxis(discMaterial->binning(), *disc);
  BOOST_CHECK_EQUAL(discAxes->getAxis(0).getNBins(), 4);
  BOOST_CHECK_EQUAL(discAxes->getAxis(0).getMin(), 10.);
  BOOST_CHECK_EQUAL(discAxes->getAxis(0).getMax(), 50.);
  BOOST_CHECK_EQUAL(discAxes->getAxis(1).getNBins(), 8);

  params.get<int>("material_binZ") = 0;
  const auto homogeneousZ = createProtoMaterial(
      params, "material", {{"binPhi", closed}, {"binZ", open}});
  BOOST_CHECK_EQUAL(homogeneousZ->binning().axisSpec(1).asEquidistant().nBins,
                    1);
}
