// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedBFieldMap_tests.cpp

// clang-format off
#define BOOST_TEST_MODULE Mapped magnetic field tests
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace tt = boost::test_tools;

using Acts::VectorHelpers::perp;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_CASE(InterpolatedBFieldMap_rz)
  {
    // definition of dummy BField
    struct BField
    {
      static Vector3D
      value(const std::array<double, 2>& rz)
      {
        double r = rz.at(0);
        double z = rz.at(1);
        // linear in r and z so interpolation should be exact
        return Vector3D(r * z, 3 * r, -2 * z);
      }
    };

    // map (x,y,z) -> (r,z)
    auto transformPos = [](const Vector3D& pos) {
      return ActsVectorD<2>(perp(pos), pos.z());
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField
        = [](const Vector3D& field, const Vector3D&) { return field; };

    // magnetic field known on grid in (r,z)
    detail::EquidistantAxis r(0.0, 4.0, 4u);
    detail::EquidistantAxis z(-5, 5, 5u);

    typedef detail::
        Grid<Vector3D, detail::EquidistantAxis, detail::EquidistantAxis>
            Grid_t;

    Grid_t g(std::make_tuple(std::move(r), std::move(z)));

    // set grid values
    for (size_t i = 1; i <= g.getNBins().at(0) + 1; ++i) {
      for (size_t j = 1; j <= g.getNBins().at(1) + 1; ++j) {
        Grid_t::index_t indices  = {{i, j}};
        const auto&     llCorner = g.getLowerLeftBinEdge(indices);
        g.at(indices)            = BField::value(llCorner);
      }
    }

    // create field mapping
    InterpolatedBFieldMap::FieldMapper<3, 2> mapper(
        transformPos, transformBField, std::move(g));
    InterpolatedBFieldMap::Config config;
    config.scale  = 1.;
    config.mapper = std::move(mapper);

    // create BField service
    InterpolatedBFieldMap b(std::move(config));

    Vector3D pos;
    pos << -3, 2.5, 1.7;
    // test the cache interface
    InterpolatedBFieldMap::Cache bCache;
    CHECK_CLOSE_REL(
        b.getField(pos, bCache), BField::value({{perp(pos), pos.z()}}), 1e-6);

    CHECK_CLOSE_REL(
        b.getField(pos, bCache), BField::value({{perp(pos), pos.z()}}), 1e-6);
    auto c = bCache.fieldCell;
    BOOST_CHECK(c.isInside(pos));
    CHECK_CLOSE_REL(
        c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

    pos << 0, 1.5, -2.5;
    InterpolatedBFieldMap::Cache bCache2;
    CHECK_CLOSE_REL(
        b.getField(pos, bCache2), BField::value({{perp(pos), pos.z()}}), 1e-6);
    c = bCache2.fieldCell;
    BOOST_CHECK(c.isInside(pos));
    CHECK_CLOSE_REL(
        c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

    pos << 2, 3, -4;
    InterpolatedBFieldMap::Cache bCache3;
    CHECK_CLOSE_REL(
        b.getField(pos, bCache3), BField::value({{perp(pos), pos.z()}}), 1e-6);
    c = bCache3.fieldCell;
    BOOST_CHECK(c.isInside(pos));
    CHECK_CLOSE_REL(
        c.getField(pos), BField::value({{perp(pos), pos.z()}}), 1e-6);

    // some field cell tests
    BOOST_CHECK(c.isInside((pos << 3, 2, -3.7).finished()));
    BOOST_CHECK(c.isInside((pos << -2, 3, -4.7).finished()));
    BOOST_CHECK(not c.isInside((pos << -2, 3, 4.7).finished()));
    BOOST_CHECK(not c.isInside((pos << 0, 2, -4.7).finished()));
    BOOST_CHECK(not c.isInside((pos << 5, 2, 14.).finished()));
  }
}  // namespace Test

}  // namespace Acts
