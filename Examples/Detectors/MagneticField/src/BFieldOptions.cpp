// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/BField/BFieldOptions.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <tuple>
#include <utility>

#include "ACTFW/Plugins/BField/BFieldUtils.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;

using InterpolatedMapper2D = Acts::InterpolatedBFieldMapper<
    Acts::detail::Grid<Acts::Vector2D, Acts::detail::EquidistantAxis,
                       Acts::detail::EquidistantAxis>>;

using InterpolatedMapper3D = Acts::InterpolatedBFieldMapper<Acts::detail::Grid<
    Acts::Vector3D, Acts::detail::EquidistantAxis,
    Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;

using InterpolatedBFieldMap2D =
    Acts::InterpolatedBFieldMap<InterpolatedMapper2D>;
using InterpolatedBFieldMap3D =
    Acts::InterpolatedBFieldMap<InterpolatedMapper3D>;

namespace FW {

namespace Options {

// common bfield options, with a bf prefix
void addBFieldOptions(boost::program_options::options_description& opt) {
  opt.add_options()("bf-map", po::value<std::string>()->default_value(""),
                    "Set this string to point to the bfield source file."
                    "That can either be a '.txt', a '.csv' or a '.root' file. "
                    "Omit for a constant magnetic field.")(
      "bf-name", po::value<std::string>()->default_value("bField"),
      "In case your field map file is given in root format, please specify "
      "the "
      "name of the TTree.")(
      "bf-gridpoints", po::value<size_t>()->default_value(100000),
      "Estimate of number of grid points, "
      "needed for allocation, only for txt and csv files.")(
      "bf-lscalor", po::value<double>()->default_value(1.),
      "The default unit for the grid "
      "points is mm. In case the grid points of your field map has another "
      "unit, please set  the scalor to mm.")(
      "bf-bscalor", po::value<double>()->default_value(1.),
      "The default unit for the magnetic field values is Tesla. In case the "
      "grid points of your field map has another unit, please set  the "
      "scalor "
      "to [T].")(
      "bf-rz", po::value<bool>()->default_value(false),
      "Please set this flag to true, if your grid points and your "
      "magnetic field values are given in 'rz'. The default is 'xyz'.")(
      "bf-foctant", po::value<bool>()->default_value(false),
      "Please set this flag to true, if your field map is only given for the "
      "first octant/quadrant and should be symmetrically created for all "
      "other "
      "octants/quadrants.")(
      "bf-values",
      po::value<read_range>()->multitoken()->default_value({0., 0., 0.}),
      "In case no magnetic field map is handed over. A constant magnetic "
      "field will be created automatically. The values can be set with this "
      "options. Please hand over the coordinates in cartesian coordinates: "
      "{Bx,By,Bz} in Tesla.")(
      "bf-context-scalable", po::value<bool>()->default_value(false),
      "This is for testing the event dependent magnetic field scaling.");
}

// create the bfield maps
BFieldVariant readBField(const boost::program_options::variables_map& vm) {
  std::string bfieldmap = "constfield";

  enum BFieldMapType { constant = 0, root = 1, text = 2 };

  std::shared_ptr<InterpolatedBFieldMap2D> map2D = nullptr;
  std::shared_ptr<InterpolatedBFieldMap3D> map3D = nullptr;
  std::shared_ptr<Acts::ConstantBField> mapConst = nullptr;
  std::shared_ptr<FW::BField::ScalableBField> mapScale = nullptr;

  int bfieldmaptype = constant;
  if (vm.count("bf-map") && vm["bf-map"].template as<std::string>() != "") {
    bfieldmap = vm["bf-map"].template as<std::string>();
    std::cout << "- read in magnetic field map: "
              << vm["bf-map"].template as<std::string>() << std::endl;
    if (bfieldmap.find(".root") != std::string::npos) {
      std::cout << "- root format for magnetic field detected" << std::endl;
      bfieldmaptype = root;
    } else if (bfieldmap.find(".txt") != std::string::npos ||
               bfieldmap.find(".csv") != std::string::npos) {
      std::cout << "- txt format for magnetic field detected" << std::endl;
      bfieldmaptype = text;
    } else {
      std::cout << "- magnetic field format could not be detected";
      std::cout << " use '.root', '.txt', or '.csv'." << std::endl;
      throw std::runtime_error("Invalid BField options");
    }
  }
  if (bfieldmaptype == text && vm.count("bf-gridpoints")) {
    std::cout << "- number of points set to: "
              << vm["bf-gridpoints"].template as<size_t>() << std::endl;
  }
  double lscalor = 1.;
  if (bfieldmaptype != constant && vm.count("bf-lscalor")) {
    lscalor = vm["bf-lscalor"].template as<double>();
    std::cout << "- length scalor to mm set to: " << lscalor << std::endl;
  }
  double bscalor = 1.;
  if (vm.count("bf-bscalor")) {
    bscalor = vm["bf-bscalor"].template as<double>();
    std::cout << "- BField (scalor to/in) Tesla set to: " << bscalor
              << std::endl;
  }
  if (bfieldmaptype != constant && vm["bf-rz"].template as<bool>())
    std::cout << "- BField map is given in 'rz' coordiantes." << std::endl;
  else if (bfieldmaptype != constant)
    std::cout << "- BField map is given in 'xyz' coordiantes." << std::endl;

  if (bfieldmaptype != constant && vm["bf-foctant"].template as<bool>()) {
    std::cout
        << "- Only the first octant/quadrant is given, bField map will be "
           "symmetrically created for all other octants/quadrants"
        << std::endl;
  }

  // Declare the mapper
  double lengthUnit = lscalor * Acts::units::_mm;
  double BFieldUnit = bscalor * Acts::units::_T;

  // set the mapper - foort
  if (bfieldmaptype == root) {
    if (vm["bf-rz"].template as<bool>()) {
      auto mapper2D = FW::BField::root::fieldMapperRZ(
          [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
            return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
          },
          vm["bf-map"].template as<std::string>(),
          vm["bf-name"].template as<std::string>(), lengthUnit, BFieldUnit,
          vm["bf-foctant"].template as<bool>());

      // create field mapping
      InterpolatedBFieldMap2D::Config config2D(std::move(mapper2D));
      config2D.scale = bscalor;
      // create BField service
      return std::make_shared<InterpolatedBFieldMap2D>(std::move(config2D));

    } else {
      auto mapper3D = FW::BField::root::fieldMapperXYZ(
          [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
            return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2)) +
                    binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
          },
          vm["bf-map"].template as<std::string>(),
          vm["bf-name"].template as<std::string>(), lengthUnit, BFieldUnit,
          vm["bf-foctant"].template as<bool>());

      // create field mapping
      InterpolatedBFieldMap3D::Config config3D(std::move(mapper3D));
      config3D.scale = bscalor;
      // create BField service
      return std::make_shared<InterpolatedBFieldMap3D>(std::move(config3D));
    }
  } else if (bfieldmaptype == text) {
    if (vm["bf-rz"].template as<bool>()) {
      auto mapper2D = FW::BField::txt::fieldMapperRZ(
          [](std::array<size_t, 2> binsRZ, std::array<size_t, 2> nBinsRZ) {
            return (binsRZ.at(1) * nBinsRZ.at(0) + binsRZ.at(0));
          },
          vm["bf-map"].template as<std::string>(), lengthUnit, BFieldUnit,
          vm["bf-gridpoints"].template as<size_t>(),
          vm["bf-foctant"].template as<bool>());

      // create field mapping
      InterpolatedBFieldMap2D::Config config2D(std::move(mapper2D));
      config2D.scale = bscalor;
      // create BField service
      return std::make_shared<InterpolatedBFieldMap2D>(std::move(config2D));

    } else {
      auto mapper3D = FW::BField::txt::fieldMapperXYZ(
          [](std::array<size_t, 3> binsXYZ, std::array<size_t, 3> nBinsXYZ) {
            return (binsXYZ.at(0) * (nBinsXYZ.at(1) * nBinsXYZ.at(2)) +
                    binsXYZ.at(1) * nBinsXYZ.at(2) + binsXYZ.at(2));
          },
          vm["bf-map"].template as<std::string>(), lengthUnit, BFieldUnit,
          vm["bf-gridpoints"].template as<size_t>(),
          vm["bf-foctant"].template as<bool>());

      // create field mapping
      InterpolatedBFieldMap3D::Config config3D(std::move(mapper3D));
      config3D.scale = bscalor;
      // create BField service
      return std::make_shared<InterpolatedBFieldMap3D>(std::move(config3D));
    }
  } else {  // constant
    // No bfield map is handed over
    // get the constant bField values
    auto bFieldValues = vm["bf-values"].template as<read_range>();
    if (bFieldValues.size() != 3) {
      throw std::invalid_argument(
          "- The values handed over for the constant magnetic field "
          "have wrong dimension. Needs to have 3 dimension. Please "
          "hand over the coordinates in cartesian coordinates: "
          "{Bx,By,Bz} in Tesla.");
    }
    if (vm["bf-context-scalable"].template as<bool>()) {
      // Create the scalable magnetic field
      return std::make_shared<FW::BField::ScalableBField>(
          bFieldValues.at(0) * Acts::units::_T,
          bFieldValues.at(1) * Acts::units::_T,
          bFieldValues.at(2) * Acts::units::_T);
    } else {
      // Create the constant magnetic field
      return std::make_shared<Acts::ConstantBField>(
          bFieldValues.at(0) * Acts::units::_T,
          bFieldValues.at(1) * Acts::units::_T,
          bFieldValues.at(2) * Acts::units::_T);
    }
  }
}
}  // namespace Options
}  // namespace FW
