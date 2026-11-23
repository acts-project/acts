// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/Detray/DetrayPropagator.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <string>

// detray includes
#include <detray/core/detector.hpp>
#include <detray/detectors/odd_metadata.hpp>
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/navigation/volume_graph.hpp>
#include <detray/propagator/line_stepper.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/propagator/rk_stepper.hpp>
#include <detray/test/common/bfield.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;
using namespace ActsPython;

namespace {

using DetrayScalar = float;
using DetrayConstField =
    typename detray::bfield::const_field_t<DetrayScalar>::view_t;

/// Bind StraightLinePropagator<tag> and HelixPropagator<tag> for a detector
template <typename detector_t>
void bindPropagators(py::module_& m, const std::string& tag) {
  using Algebra = typename detector_t::algebra_type;

  m.def(
      ("StraightLinePropagator" + tag).c_str(),
      [](std::shared_ptr<detector_t> detector,
         vecmem::memory_resource& memoryResource, bool sterile,
         Logging::Level logLevel, std::array<unsigned int, 2> searchWindow)
          -> std::shared_ptr<PropagatorInterface> {
        using Stepper = detray::line_stepper<Algebra>;
        auto logger = getDefaultLogger("DetrayPropagator", logLevel);
        if (sterile) {
          return std::make_shared<DetraySterilePropagator<Stepper, detector_t>>(
              detector, memoryResource, std::move(logger), searchWindow);
        }
        return std::make_shared<DetrayPropagator<Stepper, detector_t>>(
            detector, memoryResource, std::move(logger), searchWindow);
      },
      "detrayDetector"_a, "memoryResource"_a, "sterile"_a,
      "logLevel"_a = Logging::INFO,
      "searchWindow"_a = std::array<unsigned int, 2>{0u, 0u});

  m.def(("HelixPropagator" + tag).c_str(),
        [](std::shared_ptr<detector_t> detector,
           vecmem::memory_resource& memoryResource, bool sterile,
           std::array<double, 3> bField, Logging::Level logLevel,
           std::array<unsigned int, 2> searchWindow)
            -> std::shared_ptr<PropagatorInterface> {
          using Stepper = detray::rk_stepper<DetrayConstField, Algebra>;
          // ACTS and detray share native field units; the constant backend's
          // view holds the field value itself
          const detray::dvector3D<Algebra> b{
              static_cast<DetrayScalar>(bField[0]),
              static_cast<DetrayScalar>(bField[1]),
              static_cast<DetrayScalar>(bField[2])};
          const auto field = detray::create_const_field<DetrayScalar>(b);
          DetrayConstField view(field);
          auto logger = getDefaultLogger("DetrayPropagator", logLevel);
          if (sterile) {
            return std::make_shared<
                DetraySterilePropagator<Stepper, detector_t, DetrayConstField>>(
                detector, memoryResource, std::move(logger), searchWindow,
                view);
          }
          return std::make_shared<
              DetrayPropagator<Stepper, detector_t, DetrayConstField>>(
              detector, memoryResource, std::move(logger), searchWindow, view);
        },
        "detrayDetector"_a, "memoryResource"_a, "sterile"_a, "bField"_a,
        "logLevel"_a = Logging::INFO,
        "searchWindow"_a = std::array<unsigned int, 2>{0u, 0u});
}

}  // namespace

PYBIND11_MODULE(ActsExamplesPythonBindingsDetray, detray) {
  using DetrayDetectorODD =
      detray::host::detector<detray::odd_metadata<detray::array<float>>>;

  bindPropagators<DetrayDetectorODD>(detray, "ODD");
}
