// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryModuleHelper.hpp"

namespace dd4hep {
class Detector;
}

namespace Acts {
class TrackingGeometry;
}

namespace ActsPlugins::DD4hep::detail {
using BuildFunction = std::unique_ptr<Acts::TrackingGeometry> (*)(
    const dd4hep::Detector&, const Acts::Logger&);
const ActsGeometryModuleV1* getGeometryModule(const char* module_abi_tag,
                                              BuildFunction buildFunc);
}  // namespace ActsPlugins::DD4hep::detail

#ifdef ACTS_GEOMETRY_MODULE_ABI_TAG
#define ACTS_DEFINE_DD4HEP_GEOMETRY_MODULE(build_function) \
  ACTS_IMPL_GEOMETRY_MODULE_ENTRY(                         \
      ActsPlugins::DD4hep::detail::getGeometryModule(      \
          ACTS_GEOMETRY_MODULE_ABI_TAG, (build_function)))
#else
#define ACTS_DEFINE_DD4HEP_GEOMETRY_MODULE(build_function)           \
  static_assert(false,                                               \
                "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided via " \
                "CMake (use acts_add_dd4hep_geometry_module).")
#endif
