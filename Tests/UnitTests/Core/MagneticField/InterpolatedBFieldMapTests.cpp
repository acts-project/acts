// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedBFieldMap_tests.cpp

#include <boost/test/unit_test.hpp>

#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace tt = boost::test_tools;

using Acts::VectorHelpers::perp;

namespace Acts {

namespace Test {

// Create a test context
MagneticFieldContext mfContext = MagneticFieldContext();

BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_rz) {
  // definition of dummy BField
  struct BField {
    static Vector3 value(const std::array<double, 2>& rz) {
      double r = rz.at(0);
      double z = rz.at(1);
      // linear in r and z so interpolation should be exact
      return Vector3(r * z, 3 * r, -2 * z);
    }
  };

  // map (x,y,z) -> (r,z)
  auto transformPos = [](const Vector3& pos) {
    return Vector2(perp(pos), pos.z());
  };

  // map (Bx,By,Bz) -> (Bx,By,Bz)
  auto transformBField = [](const Vector3& field, const Vector3&) {
    return field;
  };

  // magnetic field known on grid in (r,z)
  detail::EquidistantAxis r(0.0, 4.0, 4u);
  detail::EquidistantAxis z(-5, 5, 5u);

  using Grid_t =
      detail::Grid<Vector3, detail::EquidistantAxis, detail::EquidistantAxis>;
  using Mapper_t = InterpolatedBFieldMapper<Grid_t>;
  using BField_t = InterpolatedBFieldMap<Mapper_t>;

  Grid_t g(std::make_tuple(std::move(r), std::move(z)));

  // set grid values
  for (size_t i = 1; i <= g.numLocalBins().at(0) + 1; ++i) {
    for (size_t j = 1; j <= g.numLocalBins().at(1) + 1; ++j) {
      Grid_t::index_t indices = {{i, j}};
      const auto& llCorner = g.lowerLeftBinEdge(indices);
      g.atLocalBins(indices) = BField::value(llCorner);
    }
  }

  // create field mapping
  Mapper_t mapper(transformPos, transformBField, std::move(g));
  BField_t::Config config(std::move(mapper));
  config.scale = 1.;

  // create BField service
  BField_t b(std::move(config));

  Vector3 pos;
  pos << -3, 2.5, 1.7;
  // test the cache interface
  auto bCacheAny = b.makeCache(mfContext);
  BField_t::Cache& bCache = bCacheAny.get<BField_t::Cache>();
  CHECK_CLOSE_REL(b.getField(pos, bCacheAny),
                  BField::value({{perp(pos), pos.z()}}), 1e-6);

  CHECK_CLOSE_REL(b.getField(pos, bCacheAny),
                  BField::value({{perp(pos), pos.z()}}), 1e-6);
  auto& c = *bCache.fieldCell;
  BOOST_CHECK(c.isInside(pos));
  CHECK_CLOSE_REL(c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

  pos << 0, 1.5, -2.5;
  bCacheAny = b.makeCache(mfContext);
  BField_t::Cache& bCache2 = bCacheAny.get<BField_t::Cache>();
  CHECK_CLOSE_REL(b.getField(pos, bCacheAny),
                  BField::value({{perp(pos), pos.z()}}), 1e-6);
  c = *bCache2.fieldCell;
  BOOST_CHECK(c.isInside(pos));
  CHECK_CLOSE_REL(c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

  pos << 2, 3, -4;
  bCacheAny = b.makeCache(mfContext);
  BField_t::Cache& bCache3 = bCacheAny.get<BField_t::Cache>();
  CHECK_CLOSE_REL(b.getField(pos, bCacheAny),
                  BField::value({{perp(pos), pos.z()}}), 1e-6);
  c = *bCache3.fieldCell;
  BOOST_CHECK(c.isInside(pos));
  CHECK_CLOSE_REL(c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

  // some field cell tests
  BOOST_CHECK(c.isInside((pos << 3, 2, -3.7).finished()));
  BOOST_CHECK(c.isInside((pos << -2, 3, -4.7).finished()));
  BOOST_CHECK(not c.isInside((pos << -2, 3, 4.7).finished()));
  BOOST_CHECK(not c.isInside((pos << 0, 2, -4.7).finished()));
  BOOST_CHECK(not c.isInside((pos << 5, 2, 14.).finished()));
}
}  // namespace Test

}  // namespace Acts
