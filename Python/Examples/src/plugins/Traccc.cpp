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

#include "traccc/examples/full_chain_algorithm.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/io/read_magnetic_field.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/io/read_cells.hpp"

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
                               Covfie::ConstantField::view_t>;

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
                               Covfie::InterpolatedField::view_t>;

          DetrayPropagator::Config cfg{detrayStore, sterile, std::move(ifield)};
          detrayPropagator = std::make_shared<DetrayPropagator>(cfg);

          return detrayPropagator;
        });
  }

  /// Read magnetic field from file
  {
    py::class_<traccc::magnetic_field>(traccc, "MagneticField");

    traccc.def("readMagneticField",
      [](const std::string& field_file) -> traccc::magnetic_field {
        traccc::magnetic_field field;
        traccc::io::read_magnetic_field(field, field_file);
        return field;
      },
      "field_file"_a);
  }

  /// Read detector description from files
  {
    py::class_<traccc::detector_design_description::host>(traccc, "DetectorDesignDescription");
    py::class_<traccc::detector_conditions_description::host>(traccc, "DetectorConditionsDescription");

    traccc.def("readDetectorDescription",
      [](const std::string& detector_file,
         const std::string& digitization_file,
         const std::string& conditions_file) {
        static vecmem::host_memory_resource mr;
        traccc::detector_design_description::host det_descr{mr};
        traccc::detector_conditions_description::host det_cond{mr};
        traccc::io::read_detector_description(
            det_descr, det_cond,
            detector_file, digitization_file, conditions_file,
            traccc::data_format::json);
        return std::make_pair(std::move(det_descr), std::move(det_cond));
      },
      "detector_file"_a, "digitization_file"_a, "conditions_file"_a = "");
  }

  /// Read cells from file
  {
    traccc.def("readCells",
      [](const std::string& directory,
         std::size_t event,
         const traccc::detector_conditions_description::host& det_cond) {
        static vecmem::host_memory_resource mr;
        traccc::edm::silicon_cell_collection::host cells{mr};
        static constexpr bool DEDUPLICATE = true;
        traccc::io::read_cells(
            cells, event, directory,
            traccc::getDefaultLogger("ReadCells", traccc::Logging::INFO),
            &det_cond, traccc::data_format::csv,
            DEDUPLICATE, false);
        return cells;
      },
      "directory"_a, "event"_a, "det_cond"_a);
  }

  /// Full chain algorithm bindings
  {
    using FullChain = traccc::cuda::full_chain_algorithm;

    // Input cell type
    py::class_<traccc::edm::silicon_cell_collection::host>(traccc, "SiliconCellCollection")
        .def(py::init([](){
            static vecmem::host_memory_resource mr;
            return traccc::edm::silicon_cell_collection::host{mr};
        }))
        .def("size", &traccc::edm::silicon_cell_collection::host::size);

    // Output track type
    py::class_<traccc::edm::track_collection<traccc::default_algebra>::host>(
        traccc, "TrackCollection")
        .def("size",
             [](const traccc::edm::track_collection<traccc::default_algebra>::host& t) {
               return t.size();
             });

    // The full chain algorithm
    py::class_<FullChain>(traccc, "FullChainAlgorithm")
        .def(py::init([](
            const traccc::detector_design_description::host& det_descr,
            const traccc::detector_conditions_description::host& det_cond,
            const traccc::magnetic_field& field
        ) {
            static vecmem::host_memory_resource mr;
            traccc::seedfinder_config finder_cfg;
            return std::make_unique<FullChain>(
                mr,
                traccc::clustering_config{},
                finder_cfg,
                traccc::spacepoint_grid_config{finder_cfg},
                traccc::seedfilter_config{},
                traccc::track_params_estimation_config{},
                FullChain::finding_algorithm::config_type{},
                FullChain::fitting_algorithm::config_type{},
                det_descr, det_cond, field, nullptr,
                traccc::getDefaultLogger("FullChain", traccc::Logging::INFO)
            );
        }),
        "det_descr"_a, "det_cond"_a, "field"_a)
        .def("__call__",
             [](const FullChain& alg,
                const traccc::edm::silicon_cell_collection::host& cells) {
               return alg(cells);
             },
             "cells"_a);
  }
}
