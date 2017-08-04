// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file InterpolatedBFieldMap_tests.cpp

#define BOOST_TEST_MODULE Mapped magnetic field tests
#include "ACTS/MagneticField/InterpolatedBFieldMap.hpp"
#include <boost/test/included/unit_test.hpp>
#include "ACTS/Utilities/detail/Axis.hpp"
#include "ACTS/Utilities/detail/Grid.hpp"

namespace tt = boost::test_tools;

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
      return ActsVectorD<2>(pos.perp(), pos.z());
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField
        = [](const Vector3D& field, const Vector3D&) { return field; };

    // magnetic field known on grid in (r,z)
    detail::EquidistantAxis r(0.0, 4.0, 4u);
    detail::EquidistantAxis z(-5, 5, 5u);

    typedef detail::Grid<Vector3D,
                         detail::EquidistantAxis,
                         detail::EquidistantAxis>
        Grid_t;

    Grid_t g(std::make_tuple(std::move(r), std::move(z)));

    // set grid values
    for (size_t i = 1; i <= g.getNBins<0u>() + 1; ++i) {
      for (size_t j = 1; j <= g.getNBins<1u>() + 1; ++j) {
        Grid_t::index_t indices  = {i, j};
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
    auto c = b.getFieldCell(pos);
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
    BOOST_TEST(c.isInside(pos));
    BOOST_TEST(c.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));

    pos << 0, 1.5, -2.5;
    c = b.getFieldCell(pos);
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
    BOOST_TEST(c.isInside(pos));
    BOOST_TEST(c.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));

    pos << 2, 3, -4;
    c = b.getFieldCell(pos);
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
    BOOST_TEST(c.isInside(pos));
    BOOST_TEST(c.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));

    // some field cell tests
    BOOST_TEST(c.isInside((pos << 3, 2, -3.7).finished()));
    BOOST_TEST(c.isInside((pos << -2, 3, -4.7).finished()));
    BOOST_TEST(not c.isInside((pos << -2, 3, 4.7).finished()));
    BOOST_TEST(not c.isInside((pos << 0, 2, -4.7).finished()));
    BOOST_TEST(not c.isInside((pos << 5, 2, 14.).finished()));
  }
}  // namespace Test

}  // namespace Acts
