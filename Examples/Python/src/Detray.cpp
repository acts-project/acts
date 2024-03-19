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

#include <fstream>
#include <initializer_list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
            mex.def("DetrayConverter",
                    [](const Acts::GeometryContext& gctx,
                    const Acts::Experimental::Detector& acts_detector,
                    const std::string& name) -> bool {
                        
                        //DETRAY
                        //build a mini detector
                        typename detector_t::name_map names{};
                        vecmem::host_memory_resource host_mr;
                        detector_builder<default_metadata> det_builder{};
                        detray::io::detail::detector_components_reader<detector_t> readers;

                        readers.set_detector_name(acts_detector.name());
                        //readers.read(det_builder, names);
                        const detector_t d = det_builder.build(host_mr);
                        //bool ret2 = detray_converter(d);


                        //ACTS
                        //call the converters from plugin 
                        //bool ret = detray_converter(acts_detector);
                        //const detector_t d_detray = ;                           

                        return (detray_tree_converter(acts_detector, gctx));
                    });
        }
    }
}