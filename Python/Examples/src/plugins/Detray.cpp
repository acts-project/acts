// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detray/DetrayPropagator.hpp"
#include "ActsExamples/Detray/DetrayStore.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

// detray includes
#include <detray/core/detector.hpp>
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/navigation/volume_graph.hpp>
#include <detray/propagator/line_stepper.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/propagator/rk_stepper.hpp>
#include <pybind11/pybind11.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsDetray, detray) {
  /// Define host detray store
  {
    py::class_<DetrayHostStore, std::shared_ptr<DetrayHostStore>>(
        detray, "DetrayHostStore");

    detray.def("readDetectorHost", [](const std::string& geometry,
                                      const std::string& materials,
                                      const std::string& grids) {
      auto mr = std::make_shared<vecmem::host_memory_resource>();
      auto reader_cfg = detray::io::detector_reader_config{};
      reader_cfg.add_file(geometry);
      if (!materials.empty())
        reader_cfg.add_file(materials);
      if (!grids.empty())
        reader_cfg.add_file(grids);
      auto [det, names] =
          detray::io::read_detector<DetrayHostDetector>(*mr, reader_cfg);
      return DetrayHostStore{std::move(mr), std::move(det)};
    });
  }

  /// Propagators
  {
    detray.def(
        "createSlPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           bool sterile) {
          using DetrayLineStepper =
              detray::line_stepper<typename DetrayHostDetector::algebra_type>;
          using DP = DetrayPropagator<DetrayLineStepper, DetrayHostStore>;
          DP::Config cfg{detrayStore, sterile};
          return std::shared_ptr<PropagatorInterface>(
              std::make_shared<DP>(cfg));
        },
        "store"_a, "sterile"_a = false);

    detray.def(
        "createConstBFieldPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           Covfie::ConstantField cfield, bool sterile) {
          using Stepper =
              detray::rk_stepper<Covfie::ConstantField::view_t,
                                 typename DetrayHostDetector::algebra_type>;
          using DP = DetrayPropagator<Stepper, DetrayHostStore,
                                      Covfie::ConstantField::view_t>;
          DP::Config cfg{detrayStore, sterile, cfield};
          return std::shared_ptr<PropagatorInterface>(
              std::make_shared<DP>(cfg));
        },
        "store"_a, "field"_a, "sterile"_a = false);
  }
}
