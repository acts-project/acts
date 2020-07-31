// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <limits>
#include <vector>

static constexpr auto eps = 2 * std::numeric_limits<float>::epsilon();

BOOST_AUTO_TEST_SUITE(material_properties)

BOOST_AUTO_TEST_CASE(construct_simple) {
  // construct from separate arguments
  Acts::MaterialProperties fromArgs(1., 2., 3., 4., 5., 6.);
  /// construct from material and thickness
  Acts::MaterialProperties fromMaterial(Acts::Material(1., 2., 3., 4., 5.), 6.);

  CHECK_CLOSE_REL(fromArgs.thickness(), 6., eps);
  CHECK_CLOSE_REL(fromArgs.thicknessInX0(), 6., eps);
  CHECK_CLOSE_REL(fromArgs.thicknessInL0(), 3., eps);
  CHECK_CLOSE_REL(fromMaterial.thickness(), 6., eps);
  CHECK_CLOSE_REL(fromMaterial.thicknessInX0(), 6., eps);
  CHECK_CLOSE_REL(fromMaterial.thicknessInL0(), 3., eps);
  BOOST_CHECK_EQUAL(fromArgs.material(), fromMaterial.material());
  BOOST_CHECK_EQUAL(fromArgs, fromMaterial);
}

static const Acts::MaterialProperties a(1., 2., 3., 4., 5., 1.);
static const Acts::MaterialProperties b(2., 4., 6., 8., 10., 2.);
static const Acts::MaterialProperties c(4., 8., 12., 16., 20., 3.);
static const std::vector<Acts::MaterialProperties> components = {a, b, c};

BOOST_AUTO_TEST_CASE(construct_compound) {
  Acts::MaterialProperties abc(components);

  // consistency checks
  CHECK_CLOSE_REL(abc.thickness() / abc.material().X0(), abc.thicknessInX0(),
                  eps);
  CHECK_CLOSE_REL(abc.thickness() / abc.material().L0(), abc.thicknessInL0(),
                  eps);

  // absolute and relative thicknesses are additive
  CHECK_CLOSE_REL(abc.thickness(),
                  a.thickness() + b.thickness() + c.thickness(), eps);
  CHECK_CLOSE_REL(abc.thicknessInX0(),
                  a.thicknessInX0() + b.thicknessInX0() + c.thicknessInX0(),
                  eps);
  CHECK_CLOSE_REL(abc.thicknessInL0(),
                  a.thicknessInL0() + b.thicknessInL0() + c.thicknessInL0(),
                  eps);
  // The density scales with the thickness
  CHECK_CLOSE_REL(abc.material().massDensity(),
                  (a.thickness() * a.material().massDensity() +
                   b.thickness() * b.material().massDensity() +
                   c.thickness() * c.material().massDensity()) /
                      (a.thickness() + b.thickness() + c.thickness()),
                  eps);
}

BOOST_AUTO_TEST_CASE(scale_thickness) {
  Acts::MaterialProperties mat(1., 2., 3., 4., 5., 0.1);
  Acts::MaterialProperties halfMat(1., 2., 3., 4., 5., 0.05);
  Acts::MaterialProperties halfScaled = mat;
  halfScaled.scaleThickness(0.5);

  BOOST_CHECK_NE(mat, halfMat);
  BOOST_CHECK_EQUAL(halfMat, halfScaled);
  CHECK_CLOSE_REL(mat.thicknessInX0(), 2 * halfMat.thicknessInX0(), eps);
  CHECK_CLOSE_REL(mat.thicknessInL0(), 2 * halfMat.thicknessInL0(), eps);
  CHECK_CLOSE_REL(mat.thickness() * mat.material().massDensity(),
                  2. * halfMat.thickness() * halfMat.material().massDensity(),
                  eps);
}

BOOST_AUTO_TEST_SUITE_END()
