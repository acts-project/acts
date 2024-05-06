// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "detray/builders/detector_builder.hpp"
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/io/frontend/implementation/json_readers.hpp"
#include "detray/io/frontend/utils/detector_components_reader.hpp"
#include "detray/utils/consistency_checker.hpp"

#include "Detray.hpp"

#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

#include <fstream>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

namespace Acts {
class IMaterialDecorator;
}  // namespace Acts
namespace ActsExamples {
class IMaterialWriter;
class IWriter;
}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace detray;
using namespace detray::io::detail;

using detector_t = detector<default_metadata>;


//example pybind lamda function
namespace Acts::Python {

    void addDetray(Context& ctx) {

        auto [m, mex] = ctx.get("main", "examples");

        {
            py::class_<detector<default_metadata>, std::shared_ptr<detector<default_metadata>>>(m, "detray_detector");
                //.def(py::init<detector<default_metadata>&&>(), py::arg("detector"), py::return_value_policy::move);
                //.def(py::init<vecmem::memory_resource &>(),py::arg("resource"), py::return_value_policy::move)//, py::return_value_policy::move
        }
        
        {
            mex.def("DetrayPrinter", &detray::detray_detector_print);
        }

        //detray from json function-> give good json and get detector from that
            //check if my detector building(det payload) is wrong or not
            //what consistency check is catching/missing
        
        {
            mex.def("DetrayConverter",
                    [](const Acts::GeometryContext& gctx,
                    const Acts::Experimental::Detector& acts_detector,
                    const std::string& name) -> auto {//detector_t
                        
                        //auto d_detray = std::make_shared<detector_t>(detray_from_json());
                        vecmem::host_memory_resource host_mr;
                        //auto d_detray = detray_from_json(host_mr);
                        
                        auto d_detray = detray_tree_converter(acts_detector, gctx, host_mr);   
                        //bool res = detray_tree_converter(acts_detector, gctx);
                        
                        //auto d_detray = std::make_shared<detector_t>(detray_tree_converter(acts_detector, gctx));
                        //detector_t d_detray(std::move(detray_tree_converter(acts_detector, gctx)));
                        //auto d_detray = std::make_shared<detector_t>(detray_tree_converter(acts_detector, gctx));
                        return true;
                        //return std::move(detray_tree_converter(acts_detector, gctx));
                        //return detray_tree_converter(acts_detector, gctx);
                        //return true;
                    });//, py::return_value_policy::move
        }
    }
}