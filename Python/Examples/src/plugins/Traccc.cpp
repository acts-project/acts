// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Traccc/DetrayPropagator.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

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

PYBIND11_MODULE(ActsExamplesPythonBindingsTraccc, traccc) {
  /// Define host detray store
  {
    py::class_<DetrayHostStore, std::shared_ptr<DetrayHostStore>>(
        traccc, "DetrayHostStore");

    /// Read the detray detector from files
    /// @param geometry the geometry file name
    /// @param materials the material file name
    /// @param grids the surface grids file name
    traccc.def("readDetectorHost", [](const std::string& geometry,
                                      const std::string& materials,
                                      const std::string& grids) {
      auto mr = std::make_shared<vecmem::host_memory_resource>();

      auto reader_cfg = detray::io::detector_reader_config{};
      reader_cfg.add_file(geometry);
      if (materials.empty() == false) {
        reader_cfg.add_file(materials);
      }
      if (grids.empty() == false) {
        reader_cfg.add_file(grids);
      }

      // Read the json files
      auto [det, names] =
          detray::io::read_detector<DetrayHostDetector>(*mr, reader_cfg);
      return DetrayHostStore{std::move(mr), std::move(det)};
    });
  }

  /// Define the DetrayPropagator straight line propagator
  {
    traccc.def(
        "createSlPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           bool sterile = false) {
          std::shared_ptr<PropagatorInterface> detrayPropagator = nullptr;

          using DetrayLineStepper =
              detray::line_stepper<typename DetrayHostDetector::algebra_type>;

          using DetrayPropagator =
              DetrayPropagator<DetrayLineStepper, DetrayHostStore>;

          DetrayPropagator::Config cfg{detrayStore, sterile};
          detrayPropagator = std::make_shared<DetrayPropagator>(cfg);
          return detrayPropagator;
        });
  }

  /// Define the DetrayPropagator with a covfie constant b field
  {
    traccc.def(
        "createConstBFieldPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           Covfie::ConstantField cfield, bool sterile = false) {
          std::shared_ptr<PropagatorInterface> detrayPropagator = nullptr;

          // Runge-Kutta-Nystrom stepper (field integration)
          using DetrayRknStepper =
              detray::rk_stepper<Covfie::ConstantField::view_t,
                                 typename DetrayHostDetector::algebra_type>;

          using DetrayPropagator =
              DetrayPropagator<DetrayRknStepper, DetrayHostStore,
                               Covfie::ConstantField>;

          DetrayPropagator::Config cfg{detrayStore, sterile, cfield};
          detrayPropagator = std::make_shared<DetrayPropagator>(cfg);
          return detrayPropagator;
        });
  }

  /// Define the DetrayPropagator with a covfie interpolated b field
  {
    traccc.def(
        "createInterpolatedBFieldPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           Covfie::InterpolatedField ifield, bool sterile = false) {
          std::shared_ptr<PropagatorInterface> detrayPropagator = nullptr;

          // Runge-Kutta-Nystrom stepper (field integration)
          using DetrayRknStepper =
              detray::rk_stepper<Covfie::InterpolatedField::view_t,
                                 typename DetrayHostDetector::algebra_type>;

          using DetrayPropagator =
              DetrayPropagator<DetrayRknStepper, DetrayHostStore,
                               Covfie::InterpolatedField>;

          DetrayPropagator::Config cfg{detrayStore, sterile, std::move(ifield)};
          detrayPropagator = std::make_shared<DetrayPropagator>(cfg);

          return detrayPropagator;
        });
  }
  /// Define the DetrayPropagator with interpolated b field
}
