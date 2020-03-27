// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <string>

#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "BFieldWritingBase.hpp"

namespace po = boost::program_options;

template <class T>
struct always_false : std::false_type {};

/// The main executable
///
/// Creates an InterpolatedBFieldMap from a txt or csv file and writes out the
/// grid points and values of the map into root format. The Field can then be
/// displayed using the root script printBField.cpp

int main(int argc, char* argv[]) {
  using boost::program_options::value;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addBFieldOptions(desc);
  desc.add_options()("bf-file-out",
                     value<std::string>()->default_value("BFieldOut.root"),
                     "Set this name for an output root file.")(
      "bf-map-out", value<std::string>()->default_value("bField"),
      "Set this name for the tree in the out file.")(
      "bf-out-rz", value<bool>()->default_value(false),
      "Please set this flag to true, if you want to print out the field map in "
      "cylinder coordinates (r,z). The default are cartesian coordinates "
      "(x,y,z). ")(
      "bf-rRange", value<read_range>()->multitoken(),
      "[optional] range which the bfield map should be written out in either r "
      "(cylinder "
      "coordinates) or x/y (cartesian coordinates)  in [mm]. In case no value "
      "is handed over the whole map will be written out. Please "
      "hand over by simply seperating the values by space")(
      "bf-zRange", value<read_range>()->multitoken(),
      "[optional] range which the bfield map should be written out in z in "
      "[mm].In case no value is handed over for 'bf-rRange' and 'bf-zRange the "
      "whole map will be written out. "
      "Please hand over by simply seperating the values by space")(
      "bf-rBins", value<size_t>()->default_value(200),
      "[optional] The number of bins in r. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-ZBins", value<size_t>()->default_value(300),
      "[optional] The number of bins in z. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given.")(
      "bf-PhiBins", value<size_t>()->default_value(100),
      "[optional] The number of bins in phi. This parameter only needs to be "
      "specified if 'bf-rRange' and 'bf-zRange' are given and 'bf-out-rz' is "
      "turned on.");
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  auto bFieldVar = FW::Options::readBField(vm);

  return std::visit(
      [&](auto& bField) -> int {
        using field_type =
            typename std::decay_t<decltype(bField)>::element_type;
        if constexpr (!std::is_same_v<field_type, InterpolatedBFieldMap2D> &&
                      !std::is_same_v<field_type, InterpolatedBFieldMap3D>) {
          std::cout << "Bfield map could not be read. Exiting." << std::endl;
          return EXIT_FAILURE;
        } else {
          FW::BField::writeField<field_type>(vm, bField);
          return EXIT_SUCCESS;
        }
      },
      bFieldVar);
}
