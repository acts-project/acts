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

    // map (x,yz) -> (r,z)
    auto transform = [](const Vector3D& pos) {
      return ActsVectorD<2>(pos.perp(), pos.z());
    };

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
    InterpolatedBFieldMap::FieldMapper<2> mapper(transform, std::move(g));
    InterpolatedBFieldMap::Config         c;
    c.scale  = 1.;
    c.mapper = std::move(mapper);

    // create BField service
    InterpolatedBFieldMap b(std::move(c));

    Vector3D pos;
    pos << -3, 2.5, 1.7;
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
    pos << 0, 1.5, -2.5;
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
    pos << 2, 3, -3;
    BOOST_TEST(b.getField(pos).isApprox(BField::value({pos.perp(), pos.z()})));
  }
}  // namespace Test

}  // namespace Acts
