// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE MagneticField Tests

#include <boost/test/included/unit_test.hpp>
// leave blank
#include <boost/test/data/test_case.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#include <boost/test/data/test_case.hpp>

#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

using Acts::VectorHelpers::phi;
using Acts::VectorHelpers::perp;

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace IntegrationTest {

  const double L          = 5.8 * Acts::units::_m;
  const double R          = (2.56 + 2.46) * 0.5 * 0.5 * Acts::units::_m;
  const size_t nCoils     = 1154;
  const double bMagCenter = 2. * Acts::units::_T;
  const size_t nBinsR     = 150;
  const size_t nBinsZ     = 200;

  // std::ofstream validfs("magfield_lookup.csv");
  // validfs << "x;y;z;B_x;B_y;B_z;Bm_x;Bm_y;Bm_z" << std::endl;

  InterpolatedBFieldMap
  makeFieldMap(const SolenoidBField& field)
  {

    std::ofstream ostr("solenoidmap.csv");
    ostr << "i;j;r;z;B_r;B_z" << std::endl;

    double rMin = -0.1;
    double rMax = R * 2.;
    double zMin = 2 * (-L / 2.);
    double zMax = 2 * (L / 2.);

    std::cout << "rMin = " << rMin << std::endl;
    std::cout << "rMax = " << rMax << std::endl;
    std::cout << "zMin = " << zMin << std::endl;
    std::cout << "zMax = " << zMax << std::endl;

    double stepZ = std::abs(zMax - zMin) / (nBinsZ - 1);
    double stepR = std::abs(rMax - rMin) / (nBinsR - 1);

    rMax += stepR;
    zMax += stepZ;

    // Create the axis for the grid
    Acts::detail::EquidistantAxis rAxis(rMin, rMax, nBinsR);
    Acts::detail::EquidistantAxis zAxis(zMin, zMax, nBinsZ);

    // Create the grid
    using Grid_t = Acts::detail::Grid<Acts::Vector2D,
                                      Acts::detail::EquidistantAxis,
                                      Acts::detail::EquidistantAxis>;
    Grid_t grid(std::make_tuple(std::move(rAxis), std::move(zAxis)));

    // [3] Create the transformation for the position
    // map (x,y,z) -> (r,z)
    auto transformPos = [](const Acts::Vector3D& pos) {
      return Acts::Vector2D(perp(pos), pos.z());
    };

    // [4] Create the transformation for the bfield
    // map (Br,Bz) -> (Bx,By,Bz)
    auto transformBField = [](const Acts::Vector2D& field,
                              const Acts::Vector3D& pos) {
      return Acts::Vector3D(
          field.x() * cos(phi(pos)), field.x() * sin(phi(pos)), field.y());
    };

    // iterate over all bins, set their value to the solenoid value
    // at their lower left position
    for (size_t i = 0; i <= nBinsR + 1; i++) {
      for (size_t j = 0; j <= nBinsZ + 1; j++) {
        // std::cout << "(i,j) = " << i << "," << j << std::endl;
        Grid_t::index_t index({i, j});
        if (i == 0 || j == 0 || i == nBinsR + 1 || j == nBinsZ + 1) {
          // under or overflow bin, set zero
          // std::cout << "-> under / overflow" << std::endl;
          grid.at(index) = Grid_t::value_type(0, 0);
        } else {
          // regular bin, get lower left boundary
          // std::cout << "-> regular bin" << std::endl;
          Grid_t::point_t lowerLeft = grid.getLowerLeftBinEdge(index);
          // std::cout << "lowerLeft = " << lowerLeft[0] << ", " << lowerLeft[1]
          // << std::endl;;
          // do lookup
          Vector2D B = field.getField(Vector2D(lowerLeft[0], lowerLeft[1]));
          // std::cout << "B = " << B[0] << ", " << B[1] << std::endl;
          grid.at(index) = B;

          ostr << i << ";" << j << ";" << lowerLeft[0] << ";" << lowerLeft[1];
          ostr << ";" << B[0] << ";" << B[1] << std::endl;
        }
      }
    }

    // [5] Create the mapper & BField Service
    // create field mapping
    Acts::InterpolatedBFieldMap::FieldMapper<2, 2> mapper(
        transformPos, transformBField, std::move(grid));

    Acts::InterpolatedBFieldMap::Config cfg;
    cfg.mapper = std::move(mapper);
    return Acts::InterpolatedBFieldMap(std::move(cfg));
  }

  Acts::SolenoidBField        bSolenoidField({R, L, nCoils, bMagCenter});
  Acts::InterpolatedBFieldMap bFieldMap(makeFieldMap(bSolenoidField));

  struct StreamWrapper
  {
    StreamWrapper(std::ofstream ofstr) : m_ofstr(std::move(ofstr))
    {
      m_ofstr << "x;y;z;B_x;B_y;B_z;Bm_x;Bm_y;Bm_z" << std::endl;
    }

    std::ofstream m_ofstr;
  };

  StreamWrapper valid(std::ofstream("magfield_lookup.csv"));

  // struct TestSetup {
  // TestSetup()
  //: bSolenoidField({R, L, nCoils, bMagCenter}),
  // bFieldMap(makeFieldMap(bSolenoidField))
  //{
  // std::cout << "global setup\n";
  //}

  // static SolenoidBField bSolenoidField;
  // static InterpolatedBFieldMap bFieldMap;
  //};

  // BOOST_GLOBAL_FIXTURE(TestSetup);

  const int ntests = 1000000;
  BOOST_DATA_TEST_CASE(
      constant_bfieldorward_propagation_,
      bdata::random((bdata::seed   = 1,
                     bdata::engine = std::mt19937(),
                     bdata::distribution
                     = std::uniform_real_distribution<>(1.5 * (-L / 2.),
                                                        1.5 * L / 2.)))
          ^ bdata::random((bdata::seed   = 2,
                           bdata::engine = std::mt19937(),
                           bdata::distribution
                           = std::uniform_real_distribution<>(0, R * 1.5)))
          ^ bdata::random((bdata::seed   = 3,
                           bdata::engine = std::mt19937(),
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::xrange(ntests),
      z,
      r,
      phi,
      index)
  {
    (void)index;
    // std::cout << z << " " << r << " " << phi << std::endl;
    if (index % 1000 == 0) {
      std::cout << index << std::endl;
    }

    Vector3D pos(r * std::cos(phi), r * std::sin(phi), z);
    Vector3D B  = bSolenoidField.getField(pos) / Acts::units::_T;
    Vector3D Bm = bFieldMap.getField(pos) / Acts::units::_T;

    std::ofstream& ofstr = valid.m_ofstr;
    ofstr << pos.x() << ";" << pos.y() << ";" << pos.z() << ";";
    ofstr << B.x() << ";" << B.y() << ";" << B.z() << ";";
    ofstr << Bm.x() << ";" << Bm.y() << ";" << Bm.z() << std::endl;
  }

}  // namespace IntegrationTest

}  // namespace Acts
