// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBFieldMap.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include "PropagationTestHelper.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {

/// Get the ATLAS field from :
/// http://acts.web.cern.ch/ACTS/data/AtlasField/AtlasField.tar.gz
/// to run this

namespace IntegrationTest {

// Create a mapper from the a text file
InterpolatedBFieldMap::FieldMapper<3, 3> readFieldXYZ(
    std::function<size_t(std::array<size_t, 3> binsXYZ,
                         std::array<size_t, 3> nBinsXYZ)>
        localToGlobalBin,
    std::string fieldMapFile = "Field.txt", double lengthUnit = 1.,
    double BFieldUnit = 1., size_t nPoints = 100000, bool firstOctant = false) {
  /// [1] Read in field map file
  // Grid position points in x, y and z
  std::vector<double> xPos;
  std::vector<double> yPos;
  std::vector<double> zPos;
  // components of magnetic field on grid points
  std::vector<Acts::Vector3> bField;
  // reserve estimated size
  xPos.reserve(nPoints);
  yPos.reserve(nPoints);
  zPos.reserve(nPoints);
  bField.reserve(nPoints);
  // [1] Read in file and fill values
  std::ifstream map_file(fieldMapFile.c_str(), std::ios::in);
  std::string line;
  double x = 0., y = 0., z = 0.;
  double bx = 0., by = 0., bz = 0.;
  while (std::getline(map_file, line)) {
    if (line.empty() || line[0] == '%' || line[0] == '#' ||
        line.find_first_not_of(' ') == std::string::npos)
      continue;

    std::istringstream tmp(line);
    tmp >> x >> y >> z >> bx >> by >> bz;
    xPos.push_back(x);
    yPos.push_back(y);
    zPos.push_back(z);
    bField.push_back(Acts::Vector3(bx, by, bz));
  }
  map_file.close();

  return fieldMapperXYZ(localToGlobalBin, xPos, yPos, zPos, bField, lengthUnit,
                        BFieldUnit, firstOctant);
}

// create a bfiel map from a mapper
std::shared_ptr<const InterpolatedBFieldMap> atlasBField(
    std::string fieldMapFile = "Field.txt") {
  // Declare the mapper
  Concepts ::AnyFieldLookup<> mapper;
  double lengthUnit = UnitConstants::mm;
  double BFieldUnit = UnitConstants::T;
  // read the field x,y,z from a text file
  mapper = readFieldXYZ(
      [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
        return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2)) +
                binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
      },
      fieldMapFile, lengthUnit, BFieldUnit);
  // create the config
  InterpolatedBFieldMap::Config config;
  config.scale = 1.;
  config.mapper = std::move(mapper);
  // make the interpolated field
  return std::make_shared<const InterpolatedBFieldMap>(std::move(config));
}

double Bz = 2_T;

using BFieldType = InterpolatedBFieldMap;
using EigenStepperType = EigenStepper<>;
using AtlasStepperType = AtlasStepper;
using EigenPropagatorType = Propagator<EigenStepperType>;
using AtlasPropagatorType = Propagator<AtlasStepperType>;

auto bField = atlasBField("Field.txt");

EigenStepperType estepper(bField);
EigenPropagatorType epropagator(std::move(estepper));
AtlasStepperType astepper(bField);
AtlasPropagatorType apropagator(std::move(astepper));

// The actual test - needs to be included to avoid
// template inside template definition through boost
#include "PropagationTestBase.hpp"

}  // namespace IntegrationTest

}  // namespace Acts
