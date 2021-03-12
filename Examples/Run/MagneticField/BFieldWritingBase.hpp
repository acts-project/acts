// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/Root/RootBFieldWriter.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace BField {

template <typename bfield_t>
void writeField(boost::program_options::variables_map vm,
                const bfield_t& bField) {
  using Writer = ActsExamples::RootBFieldWriter<bfield_t>;
  using Config = typename Writer::Config;
  using GridType = typename Writer::GridType;

  // Write the interpolated magnetic field
  Config writerConfig;
  if (vm["bf-out-rz"].template as<bool>()) {
    writerConfig.gridType = GridType::rz;
  } else {
    writerConfig.gridType = GridType::xyz;
  }
  writerConfig.treeName = vm["bf-map-out"].template as<std::string>();
  writerConfig.fileName = vm["bf-file-out"].template as<std::string>();
  writerConfig.bField = &bField;
  std::cout << "setting rBounds" << std::endl;
  if (vm.count("bf-rRange") && vm.count("bf-zRange")) {
    auto rBounds = vm["bf-rRange"].template as<Options::Reals<2>>();
    auto zBounds = vm["bf-zRange"].template as<Options::Reals<2>>();
    writerConfig.rBounds = {{rBounds[0] * Acts::UnitConstants::mm,
                             rBounds[1] * Acts::UnitConstants::mm}};
    writerConfig.zBounds = {{zBounds[0] * Acts::UnitConstants::mm,
                             zBounds[1] * Acts::UnitConstants::mm}};
  }
  writerConfig.rBins = vm["bf-rBins"].template as<size_t>();
  writerConfig.zBins = vm["bf-ZBins"].template as<size_t>();
  writerConfig.phiBins = vm["bf-PhiBins"].template as<size_t>();

  ActsExamples::RootBFieldWriter<bfield_t>::run(writerConfig);
}
}  // namespace BField
}  // namespace ActsExamples
