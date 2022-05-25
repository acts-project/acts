// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/MagneticField/FieldMapRootIo.hpp"
#include "ActsExamples/MagneticField/FieldMapTextIo.hpp"
#include "ActsExamples/MagneticField/ScalableBFieldService.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>
#include <stdexcept>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

void ActsExamples::Options::addMagneticFieldOptions(Description& desc) {
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  // avoid adding the options twice
  if (desc.find_nothrow("bf-constant-tesla", true) != nullptr) {
    return;
  }

  auto opt = desc.add_options();
  opt("bf-constant-tesla", value<Reals<3>>(),
      "Set a constant magnetic field vector in Tesla. If given, this takes "
      "preference over all other options.");
  opt("bf-scalable", bool_switch(),
      "If given, the constant field strength will be scaled differently in "
      "every event. This is for testing only.");
  opt("bf-scalable-scalor", value<double>()->default_value(1.25),
      "Scaling factor for the event-dependent field strength scaling. A unit "
      "value means that the field strength stays the same for every event.");
  opt("bf-map-file", value<std::string>(),
      "Read a magnetic field map from the given file. ROOT and text file "
      "formats are supported. Only used if no constant field is given.");
  opt("bf-map-tree", value<std::string>()->default_value("bField"),
      "Name of the TTree in the ROOT file. Only used if the field map is read "
      "from a ROOT file.");
  opt("bf-map-type", value<std::string>()->default_value("xyz"),
      "Either 'xyz' or 'rz' to define the type of the field map.");
  opt("bf-map-octantonly", bool_switch(),
      "If given, the field map is assumed to describe only the first "
      "octant/quadrant and the field is symmetrically extended to the full "
      "space.");
  opt("bf-map-lengthscale-mm", value<double>()->default_value(1.),
      "Optional length scale modifier for the field map grid. This options "
      "only needs to be set if the length unit in the field map file is not "
      "`mm`. The value must scale from the stored unit to the equivalent value "
      "in `mm`.");
  opt("bf-map-fieldscale-tesla", value<double>()->default_value(1.),
      "Optional field value scale modifier for the field map value. This "
      "option only needs to be set if the field value unit in the field map "
      "file is not `Tesla`. The value must scale from the stored unit to the "
      "equvalent value in `Tesla`.");
  opt("bf-solenoid-mag-tesla", value<double>()->default_value(0.),
      "The magnitude of a solenoid magnetic field in the center in `Tesla`. "
      "Only used "
      "if neither constant field nor a magnetic field map is given.");
  opt("bf-solenoid-length", value<double>()->default_value(6000),
      "The length of the solenoid magnetic field in `mm`.");
  opt("bf-solenoid-radius", value<double>()->default_value(1200),
      "The radius of the solenoid magnetic field in `mm`.");
  opt("bf-solenoid-ncoils", value<size_t>()->default_value(1194),
      "Number of coils for the solenoid magnetic field.");
  opt("bf-solenoid-map-rlim",
      value<Interval>()->value_name("MIN:MAX")->default_value({0, 1200}),
      "The length bounds of the grid created from the analytical solenoid "
      "field in `mm`.");
  opt("bf-solenoid-map-zlim",
      value<Interval>()->value_name("MIN:MAX")->default_value({-3000, 3000}),
      "The radius bounds of the grid created from the analytical solenoid "
      "field in `mm`.");
  opt("bf-solenoid-map-nbins", value<Reals<2>>()->default_value({{150, 200}}),
      "The number of bins in r-z directions for the grid created from the "
      "analytical solenoid field.");
}

void ActsExamples::Options::setupMagneticFieldServices(const Variables& vars,
                                                       Sequencer& sequencer) {
  if (vars["bf-scalable"].as<bool>()) {
    ScalableBFieldService::Config sbfCfg;
    sbfCfg.scalor = vars["bf-scalable-scalor"].as<double>();
    sequencer.addContextDecorator(
        std::make_shared<ScalableBFieldService>(sbfCfg, Acts::Logging::INFO));
  }
}

std::shared_ptr<Acts::MagneticFieldProvider>
ActsExamples::Options::readMagneticField(const Variables& vars) {
  using namespace ActsExamples::detail;
  using boost::filesystem::path;

  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("MagneticField", Acts::Logging::INFO));

  // first option: create a constant field
  if (vars.count("bf-constant-tesla")) {
    const auto values = vars["bf-constant-tesla"].as<Reals<3>>();
    Acts::Vector3 field(values[0] * Acts::UnitConstants::T,
                        values[1] * Acts::UnitConstants::T,
                        values[2] * Acts::UnitConstants::T);
    if (vars["bf-scalable"].as<bool>()) {
      ACTS_INFO("Use a constant magnetic field with per-event scaling");
      return std::make_shared<ScalableBField>(field);
    } else {
      ACTS_INFO("Use a constant magnetic field");
      return std::make_shared<Acts::ConstantBField>(field);
    }
  }

  // second option: read a field map from a file
  if (vars.count("bf-map-file")) {
    const path file = vars["bf-map-file"].as<std::string>();
    const auto tree = vars["bf-map-tree"].as<std::string>();
    const auto type = vars["bf-map-type"].as<std::string>();
    const auto useOctantOnly = vars["bf-map-octantonly"].as<bool>();
    const auto lengthUnit =
        vars["bf-map-lengthscale-mm"].as<double>() * Acts::UnitConstants::mm;
    const auto fieldUnit =
        vars["bf-map-fieldscale-tesla"].as<double>() * Acts::UnitConstants::T;

    bool readRoot = false;
    if (file.extension() == ".root") {
      ACTS_INFO("Read magnetic field map from ROOT file '" << file << "'");
      readRoot = true;
    } else if (file.extension() == ".txt") {
      ACTS_INFO("Read magnetic field map from text file '" << file << "'");
      readRoot = false;
    } else {
      ACTS_ERROR("'" << file
                     << "' is an unsupported magnetic field map file type");
      throw std::runtime_error("Unsupported magnetic field map file type");
    }

    if (type == "xyz") {
      auto mapBins = [](std::array<size_t, 3> bins,
                        std::array<size_t, 3> sizes) {
        return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] + bins[2]);
      };

      ACTS_INFO("Use XYZ field map");
      if (readRoot) {
        auto map = makeMagneticFieldMapXyzFromRoot(
            std::move(mapBins), file.native(), tree, lengthUnit, fieldUnit,
            useOctantOnly);
        return std::make_shared<InterpolatedMagneticField3>(std::move(map));

      } else {
        auto map = makeMagneticFieldMapXyzFromText(std::move(mapBins),
                                                   file.native(), lengthUnit,
                                                   fieldUnit, useOctantOnly);
        return std::make_shared<InterpolatedMagneticField3>(std::move(map));
      }

    } else if (type == "rz") {
      auto mapBins = [](std::array<size_t, 2> bins,
                        std::array<size_t, 2> sizes) {
        return (bins[1] * sizes[0] + bins[0]);
      };

      ACTS_INFO("Use RZ field map");
      if (readRoot) {
        auto map = makeMagneticFieldMapRzFromRoot(
            std::move(mapBins), file.native(), tree, lengthUnit, fieldUnit,
            useOctantOnly);
        return std::make_shared<InterpolatedMagneticField2>(std::move(map));

      } else {
        auto map = makeMagneticFieldMapRzFromText(std::move(mapBins),
                                                  file.native(), lengthUnit,
                                                  fieldUnit, useOctantOnly);
        return std::make_shared<InterpolatedMagneticField2>(std::move(map));
      }

    } else {
      ACTS_ERROR("'" << type << "' is an unknown magnetic field map type");
      throw std::runtime_error("Unknown magnetic field map type");
    }
  }

  // third option: create a solenoid field
  if (vars["bf-solenoid-mag-tesla"].as<double>() > 0) {
    // Construct a solenoid field
    Acts::SolenoidBField::Config solenoidConfig;
    solenoidConfig.length =
        vars["bf-solenoid-length"].as<double>() * Acts::UnitConstants::mm;
    solenoidConfig.radius =
        vars["bf-solenoid-radius"].as<double>() * Acts::UnitConstants::mm;
    solenoidConfig.nCoils = vars["bf-solenoid-ncoils"].as<size_t>();
    solenoidConfig.bMagCenter =
        vars["bf-solenoid-mag-tesla"].as<double>() * Acts::UnitConstants::T;
    ACTS_INFO("Use solenoid magnetic field with magnitude "
              << solenoidConfig.bMagCenter / Acts::UnitConstants::T
              << " Tesla at the center.");
    const auto solenoidField = Acts::SolenoidBField(solenoidConfig);
    // The parameters for creating a field map
    auto getRange = [&](const char* name, auto unit, auto& lower, auto& upper) {
      auto interval = vars[name].as<Options::Interval>();
      lower = interval.lower.value() * unit;
      upper = interval.upper.value() * unit;
    };
    std::pair<double, double> rlim, zlim;
    getRange("bf-solenoid-map-rlim", Acts::UnitConstants::mm, rlim.first,
             rlim.second);
    getRange("bf-solenoid-map-zlim", Acts::UnitConstants::mm, zlim.first,
             zlim.second);
    const auto nbins = vars["bf-solenoid-map-nbins"].as<Reals<2>>();
    auto map =
        Acts::solenoidFieldMap(rlim, zlim, {nbins[0], nbins[1]}, solenoidField);
    return std::make_shared<InterpolatedMagneticField2>(std::move(map));
  }

  // default option: no field
  ACTS_INFO("Use no magnetic field");
  return std::make_shared<Acts::NullBField>();
}
