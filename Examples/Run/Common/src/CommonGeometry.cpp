// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Geometry/MaterialWiper.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"
#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Plugins/Json/JsonMaterialDecorator.hpp>
#include <Acts/Plugins/Json/MaterialMapJsonConverter.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <string>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Geometry {

/// @brief helper method to setup the geometry
///
/// @tparam options_map_t Type of the options to be read
/// @tparam geometry_setupt_t Type of the callable geometry setup
///
/// @param vm the parsed options map
/// @param geometrySetup the callable geometry setup
///
/// @return a pair of TrackingGeometry and context decorators
std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
          std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>>
build(const boost::program_options::variables_map& vm,
      IBaseDetector& detector) {
  // Material decoration
  std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
  auto matType = vm["mat-input-type"].template as<std::string>();
  if (matType == "none") {
    matDeco = std::make_shared<const Acts::MaterialWiper>();
  } else if (matType == "file") {
    // Retrieve the filename
    auto fileName = vm["mat-input-file"].template as<std::string>();
    // json or root based decorator
    if (fileName.find(".json") != std::string::npos ||
        fileName.find(".cbor") != std::string::npos) {
      // Set up the converter first
      Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
      // Set up the json-based decorator
      matDeco = std::make_shared<const Acts::JsonMaterialDecorator>(
          jsonGeoConvConfig, fileName, Acts::Logging::INFO);
    } else if (fileName.find(".root") != std::string::npos) {
      // Set up the root-based decorator
      ActsExamples::RootMaterialDecorator::Config rootMatDecConfig;
      rootMatDecConfig.fileName = fileName;
      matDeco = std::make_shared<const ActsExamples::RootMaterialDecorator>(
          rootMatDecConfig, Acts::Logging::INFO);
    }
  }

  /// Return the geometry and context decorators
  return detector.finalize(vm, matDeco);
}

}  // namespace Geometry
}  // namespace ActsExamples
