// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Plugins/Detray/DetrayConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace detray;
using namespace detray::io::detail;

namespace Acts::Python {

void addDetray(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<detector<default_metadata>,
               std::shared_ptr<detector<default_metadata>>>(m,
                                                            "detray_detector");
  }

  { mex.def("writeToJson", &DetrayConverter::writeToJson); }
  /**
          {
              /// @brief Converts an Acts::Detector to a detray::detector
              mex.def("convertDetectorToDetray",
                      [](const Acts::GeometryContext& gctx,
                      const Acts::Experimental::Detector& acts_detector,
                      const std::string& name) -> auto
     {//detector<default_metadata>

                          // Create a host memory resource
                          vecmem::host_memory_resource host_mr;
                          // Convert Acts detector to detray detector using the
     detray converter function auto det_tuple =
     DetrayConverter::detrayConvert(acts_detector, gctx, host_mr);

                          return true; //TO DO:: cannot return tuple

                      });
          }

          **/
}
}  // namespace Acts::Python
