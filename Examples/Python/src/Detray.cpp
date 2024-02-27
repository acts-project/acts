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

//example pybind lamda function
namespace Acts::Python {

    void addDetray(Context& ctx) {

        auto [m, mex] = ctx.get("main", "examples");
        {
            mex.def("DetrayConverter",
                    [](const Acts::GeometryContext& gctx,
                    const Acts::Experimental::Detector& detector,
                    const std::string& name) -> bool {
                        //pass a boolean for checking 
                        return true;
                    });
        }
    }
}