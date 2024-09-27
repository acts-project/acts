// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Detray/DetrayConversionUtils.hpp"
#include "Acts/Plugins/Detray/DetrayConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Traccc/DetrayPropagator.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"
#include "ActsExamples/Traccc/TracccChainConfig.hpp"

#include <detray/propagator/line_stepper.hpp>
#include <pybind11/pybind11.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {

void addTracccChainHost(Context &ctx);

void addTraccc(Context &ctx) {
  auto mex = ctx.get("examples");
  auto m = mex.def_submodule("traccc");
  ctx.modules["traccc"] = m;

  /// Define host detray store
  {
    py::class_<DetrayHostStore, std::shared_ptr<DetrayHostStore>>(
        m, "DetrayHostStore");

    /// Convert the detector and create a DetrayHostStore
    ///
    /// @param gctx the geometry context
    /// @param detector the detector to be converted
    /// @param options the conversion options
    m.def("convertDetectorHost", [](const GeometryContext &gctx,
                                    const Experimental::Detector &detector,
                                    DetrayConverter::Options options) {
      return DetrayHostStore::create(gctx, detector, options);
    });
  }

  /// Define the DetrayPropagator
  {
    m.def("createSlPropagatorHost",
          [](std::shared_ptr<const DetrayHostStore> detrayStore,
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

  {
    using Config = traccc::seedfinder_config;

    auto c = py::class_<Config>(m, "SeedFinderConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(rMin);
    ACTS_PYTHON_MEMBER(collisionRegionMin);
    ACTS_PYTHON_MEMBER(collisionRegionMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(sigmaScattering);
    ACTS_PYTHON_MEMBER(maxPtScattering);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_MEMBER(beamPos);
    ACTS_PYTHON_MEMBER(radLengthPerSeed);
    ACTS_PYTHON_MEMBER(zAlign);
    ACTS_PYTHON_MEMBER(rAlign);
    ACTS_PYTHON_MEMBER(sigmaError);
    ACTS_PYTHON_MEMBER(highland);
    ACTS_PYTHON_MEMBER(maxScatteringAngle2);
    ACTS_PYTHON_MEMBER(pTPerHelixRadius);
    ACTS_PYTHON_MEMBER(minHelixDiameter2);
    ACTS_PYTHON_MEMBER(pT2perRadius);
    ACTS_PYTHON_MEMBER(phiBinDeflectionCoverage);
    ACTS_PYTHON_MEMBER(neighbor_scope);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = traccc::spacepoint_grid_config;

    auto c = py::class_<Config>(m, "SpacePointGridConfig")
                 .def(py::init<const traccc::seedfinder_config &>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(bFieldInZ);
    ACTS_PYTHON_MEMBER(minPt);
    ACTS_PYTHON_MEMBER(rMax);
    ACTS_PYTHON_MEMBER(zMax);
    ACTS_PYTHON_MEMBER(zMin);
    ACTS_PYTHON_MEMBER(deltaRMax);
    ACTS_PYTHON_MEMBER(cotThetaMax);
    ACTS_PYTHON_MEMBER(impactMax);
    ACTS_PYTHON_MEMBER(phiMin);
    ACTS_PYTHON_MEMBER(phiMax);
    ACTS_PYTHON_MEMBER(phiBinDeflectionCoverage);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = traccc::seedfilter_config;

    auto c = py::class_<Config>(m, "SeedFilterConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(deltaInvHelixDiameter);
    ACTS_PYTHON_MEMBER(impactWeightFactor);
    ACTS_PYTHON_MEMBER(compatSeedWeight);
    ACTS_PYTHON_MEMBER(deltaRMin);
    ACTS_PYTHON_MEMBER(maxSeedsPerSpM);
    ACTS_PYTHON_MEMBER(compatSeedLimit);
    ACTS_PYTHON_MEMBER(max_triplets_per_spM);
    ACTS_PYTHON_MEMBER(good_spB_min_radius);
    ACTS_PYTHON_MEMBER(good_spB_weight_increase);
    ACTS_PYTHON_MEMBER(good_spT_max_radius);
    ACTS_PYTHON_MEMBER(good_spT_weight_increase);
    ACTS_PYTHON_MEMBER(good_spB_min_weight);
    ACTS_PYTHON_MEMBER(seed_min_weight);
    ACTS_PYTHON_MEMBER(spB_min_radius);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = traccc::finding_config;

    auto c = py::class_<Config>(m, "FindingConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(max_num_branches_per_seed);
    ACTS_PYTHON_MEMBER(max_num_branches_per_surface);
    ACTS_PYTHON_MEMBER(min_track_candidates_per_track);
    ACTS_PYTHON_MEMBER(max_track_candidates_per_track);
    ACTS_PYTHON_MEMBER(max_num_skipping_per_cand);
    ACTS_PYTHON_MEMBER(min_step_length_for_next_surface);
    ACTS_PYTHON_MEMBER(max_step_counts_for_next_surface);
    ACTS_PYTHON_MEMBER(chi2_max);
    ACTS_PYTHON_MEMBER(propagation);
    ACTS_PYTHON_MEMBER(n_measurements_per_thread);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = traccc::fitting_config;

    auto c = py::class_<Config>(m, "FittingConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(n_iterations);
    ACTS_PYTHON_MEMBER(propagation);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using Config = traccc::greedy_ambiguity_resolution_algorithm::config_t;

    auto c = py::class_<Config>(m, "GreedyAmbiguityResolutionAlgorithmConfig")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(maximum_shared_hits);
    ACTS_PYTHON_MEMBER(maximum_iterations);
    ACTS_PYTHON_MEMBER(n_measurements_min);
    ACTS_PYTHON_MEMBER(check_obvious_errs);
    ACTS_PYTHON_MEMBER(measurement_id_0_warning_threshold);
    ACTS_PYTHON_MEMBER(verbose_error);
    ACTS_PYTHON_MEMBER(verbose_warning);
    ACTS_PYTHON_MEMBER(verbose_info);
    ACTS_PYTHON_MEMBER(verbose_debug);
    ACTS_PYTHON_STRUCT_END();
  }

  addTracccChainHost(ctx);
}
}  // namespace Acts::Python
