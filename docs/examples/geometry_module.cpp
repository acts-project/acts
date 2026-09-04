// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleLoader.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <filesystem>

#ifdef ACTS_BUILD_PLUGIN_DD4HEP
#include "ActsPlugins/DD4hep/GeometryModuleLoader.hpp"

#include <DD4hep/Detector.h>
#endif

void exampleLoadPlainModule(const Acts::Logger& logger,
                            const std::filesystem::path& modulePath) {
  //! [Load Plain Module]
  auto geometry = Acts::loadGeometryModule(modulePath, logger);
  // 'geometry' is a std::shared_ptr<Acts::TrackingGeometry>.
  // The shared library stays loaded until 'geometry' is destroyed.
  //! [Load Plain Module]
  (void)geometry;
}

#ifdef ACTS_BUILD_PLUGIN_DD4HEP
void exampleLoadDD4hepModule(const dd4hep::Detector& detector,
                             const Acts::Logger& logger,
                             const std::filesystem::path& modulePath) {
  //! [Load DD4hep Module]
  auto geometry = Acts::loadDD4hepGeometryModule(modulePath, detector, logger);
  // 'geometry' is a std::shared_ptr<Acts::TrackingGeometry>.
  // The shared library stays loaded until 'geometry' is destroyed.
  //! [Load DD4hep Module]
  (void)geometry;
}
#endif
