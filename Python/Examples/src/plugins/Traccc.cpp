// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccBackend.hpp"
#include "ActsExamples/Traccc/TracccCopyAlgs.hpp"
#include "ActsExamples/Traccc/TracccCopyTool.hpp"
#include "ActsExamples/Traccc/TracccDetector.hpp"
#include "ActsExamples/Traccc/TracccDeviceMeasurementSortingAlg.hpp"
#include "ActsExamples/Traccc/TracccDeviceSpacepointFormationAlg.hpp"
#include "ActsExamples/Traccc/TracccDeviceTrackFindingAlg.hpp"
#include "ActsExamples/Traccc/TracccDeviceTrackParamEstimationAlg.hpp"
#include "ActsExamples/Traccc/TracccDeviceTripletSeedingAlg.hpp"
#include "ActsExamples/Traccc/TracccEventData.hpp"
#include "ActsExamples/Traccc/TracccMagneticField.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.hpp"
#include "ActsPlugins/Detray/DetraySurfaceIdMap.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#ifdef ACTS_EXAMPLES_TRACCC_CUDA
#include "ActsExamples/Traccc/Cuda/CudaTracccBackend.hpp"
#endif

#include <array>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>
#include <traccc/bfield/magnetic_field.hpp>
#include <traccc/geometry/detector.hpp>
#include <traccc/geometry/host_detector.hpp>
#include <traccc/io/read_magnetic_field.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;
using namespace ActsPython;

namespace {

void addConfigs(py::module_& m) {
  {
    using C = traccc::seedfinder_config;
    auto c = py::class_<C>(m, "SeedFinderConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, zMin, zMax, rMin, rMax, collisionRegionMin,
                       collisionRegionMax, minPt, cotThetaMax, deltaRMin,
                       deltaRMax, deltaZMax, impactMax, sigmaScattering,
                       maxPtScattering, radLengthPerSeed, bFieldInZ);
  }
  {
    using C = traccc::seedfilter_config;
    auto c = py::class_<C>(m, "SeedFilterConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, deltaInvHelixDiameter, impactWeightFactor,
                       compatSeedWeight, deltaRMin, compatSeedLimit,
                       seed_min_weight, spB_min_radius);
  }
  {
    using C = traccc::track_params_estimation_config;
    auto c = py::class_<C>(m, "TrackParamsEstimationConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, initial_sigma, initial_sigma_qopt,
                       initial_sigma_pt_rel, initial_inflation);
  }
  {
    using C = traccc::finding_config;
    auto c = py::class_<C>(m, "FindingConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, max_num_branches_per_seed, max_num_branches_per_surface,
        min_track_candidates_per_track, max_track_candidates_per_track,
        max_num_skipping_per_cand, max_num_consecutive_skipped,
        max_num_tracks_per_measurement, min_measurement_voting_fraction,
        min_step_length_for_next_surface, max_step_counts_for_next_surface,
        chi2_max, min_p, min_pT, duplicate_removal_minimum_length,
        initial_links_per_seed);
    c.def_property(
        "search_window",
        [](const C& self) {
          return std::array<unsigned int, 2>{
              self.propagation.navigation.search_window[0],
              self.propagation.navigation.search_window[1]};
        },
        [](C& self, const std::array<unsigned int, 2>& window) {
          self.propagation.navigation.search_window = {window[0], window[1]};
        });
    c.def_property(
        "overstep_tolerance",
        [](const C& self) {
          return self.propagation.navigation.intersection.overstep_tolerance;
        },
        [](C& self, float value) {
          self.propagation.navigation.intersection.overstep_tolerance = value;
        });
  }
}

void addBackends(py::module_& m) {
  py::class_<TracccBackend, std::shared_ptr<TracccBackend>>(m, "TracccBackend")
      .def_property_readonly("name", &TracccBackend::name);

#ifdef ACTS_EXAMPLES_TRACCC_CUDA
  py::class_<CudaTracccBackend, TracccBackend,
             std::shared_ptr<CudaTracccBackend>>(m, "CudaTracccBackend")
      .def(py::init([](int device, Acts::Logging::Level level) {
             CudaTracccBackend::Config cfg;
             cfg.device = device;
             return std::make_shared<CudaTracccBackend>(
                 cfg, Acts::getDefaultLogger("CudaTracccBackend", level));
           }),
           "device"_a = 0, "level"_a = Acts::Logging::INFO);
#endif
}

void addDetector(py::module_& m) {
  py::class_<traccc::host_detector, std::shared_ptr<traccc::host_detector>>(
      m, "TracccHostDetector");

  // Convert the ACTS tracking geometry of the ODD into a type erased detray
  // detector of the traccc ODD type, returned together with the ACTS <->
  // detray surface identifier map of the conversion. The detray detector is
  // move only, so it is moved into the type erased container right after the
  // conversion. The detector keeps the memory resource it was created with
  // alive.
  m.def(
      "convertToTracccOddDetector",
      [](const ActsPlugins::DetrayGeometryConverter& converter,
         const std::shared_ptr<vecmem::memory_resource>& mr,
         const Acts::GeometryContext& gctx,
         const std::shared_ptr<const Acts::TrackingGeometry>& trackingGeometry) {
        if (mr == nullptr) {
          throw std::invalid_argument(
              "convertToTracccOddDetector: memory resource is not set");
        }
        auto geometry = converter.convert<traccc::odd_detector::metadata>(
            *mr, gctx, trackingGeometry);
        std::shared_ptr<traccc::host_detector> hostDetector(
            new traccc::host_detector(),
            [mr](traccc::host_detector* detector) { delete detector; });
        hostDetector->set<traccc::odd_detector>(std::move(*geometry.detector));
        auto surfaceIds = std::make_shared<ActsPlugins::DetraySurfaceIdMap>(
            std::move(geometry.surfaceIds));
        return std::make_pair(hostDetector, surfaceIds);
      },
      "converter"_a, "mr"_a, "gctx"_a, "trackingGeometry"_a);

  auto d = py::class_<TracccDetector, std::shared_ptr<TracccDetector>>(
               m, "TracccDetector")
               .def(py::init([](const TracccDetector::Config& cfg,
                                Acts::Logging::Level level) {
                      return std::make_shared<TracccDetector>(
                          cfg, Acts::getDefaultLogger("TracccDetector", level));
                    }),
                    "config"_a, "level"_a = Acts::Logging::INFO)
               .def_property_readonly("config", &TracccDetector::config);

  auto c = py::class_<TracccDetector::Config>(d, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, backend, trackingGeometry, detrayDetector, surfaceIds);
}

void addMagneticField(py::module_& m) {
  py::class_<traccc::magnetic_field, std::shared_ptr<traccc::magnetic_field>>(
      m, "TracccHostMagneticField");

  // Build the traccc host field from an existing covfie field, or read it
  // from a traccc (covfie) field file
  m.def(
      "makeTracccHostMagneticField",
      [](const ActsPlugins::Covfie::ConstantField& field) {
        return std::make_shared<traccc::magnetic_field>(
            ActsPlugins::Covfie::ConstantField(field));
      },
      "field"_a);
  m.def(
      "readTracccHostMagneticField",
      [](const std::filesystem::path& path) {
        auto field = std::make_shared<traccc::magnetic_field>();
        traccc::io::read_magnetic_field(*field, path.string());
        return field;
      },
      "path"_a);

  auto f =
      py::class_<TracccMagneticField, std::shared_ptr<TracccMagneticField>>(
          m, "TracccMagneticField")
          .def(py::init([](const TracccMagneticField::Config& cfg,
                           Acts::Logging::Level level) {
                 return std::make_shared<TracccMagneticField>(
                     cfg, Acts::getDefaultLogger("TracccMagneticField", level));
               }),
               "config"_a, "level"_a = Acts::Logging::INFO)
          .def_property_readonly("config", &TracccMagneticField::config);

  auto c =
      py::class_<TracccMagneticField::Config>(f, "Config").def(py::init<>());
  ACTS_PYTHON_STRUCT(c, backend, field);
}

void addCopyTool(py::module_& m) {
  py::class_<TracccCopyTool, std::shared_ptr<TracccCopyTool>>(m,
                                                              "TracccCopyTool")
      .def(py::init<std::shared_ptr<const TracccBackend>>(), "backend"_a);
}

void addAlgorithms(py::module_& m) {
  // Host <-> device copies, one instance per traccc collection
  {
    using Alg = TracccHostToDeviceCopyAlg<traccc::edm::measurement_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccHostToDeviceMeasurementCopyAlg",
                                  copyTool, input, output);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::measurement_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccDeviceToHostMeasurementCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::spacepoint_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostSpacepointCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<traccc::edm::seed_collection>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostSeedCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<
        traccc::bound_track_parameters_collection_types>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m,
                                  "TracccDeviceToHostTrackParametersCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }
  {
    using Alg = TracccDeviceToHostCopyAlg<TracccTrackContainer>;
    ACTS_PYTHON_DECLARE_ALGORITHM(Alg, m, "TracccDeviceToHostTrackCopyAlg",
                                  copyTool, input, output, inputMeasurements);
  }

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccDeviceMeasurementSortingAlg, m,
                                "TracccDeviceMeasurementSortingAlg", detector,
                                inputMeasurements, outputMeasurements);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccDeviceSpacepointFormationAlg, m,
                                "TracccDeviceSpacepointFormationAlg", detector,
                                inputMeasurements, outputSpacepoints);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccDeviceTripletSeedingAlg, m,
                                "TracccDeviceTripletSeedingAlg", detector,
                                inputSpacepoints, outputSeeds, seedFinder,
                                seedFilter);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccDeviceTrackParamEstimationAlg, m,
                                "TracccDeviceTrackParamEstimationAlg", detector,
                                magneticField, inputMeasurements, inputSpacepoints,
                                inputSeeds, outputTrackParameters, estimation);

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccDeviceTrackFindingAlg, m,
                                "TracccDeviceTrackFindingAlg", detector,
                                magneticField, inputMeasurements, inputTrackParameters,
                                outputTracks, finding);

}

}  // namespace

PYBIND11_MODULE(ActsExamplesPythonBindingsTraccc, traccc) {
  addConfigs(traccc);
  addBackends(traccc);
  addDetector(traccc);
  addMagneticField(traccc);
  addCopyTool(traccc);
  addAlgorithms(traccc);
}
