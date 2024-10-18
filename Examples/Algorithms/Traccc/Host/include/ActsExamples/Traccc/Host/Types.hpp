// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/field.hpp>

#include "detray/core/detector.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

namespace ActsExamples::Traccc::Host {

struct Types {
  using FieldType =
      covfie::field<covfie::backend::constant<covfie::vector::float3,
                                              covfie::vector::float3>>;
  using DetectorType =
      detray::detector<detray::default_metadata, detray::host_container_types>;
  using StepperType =
      detray::rk_stepper<FieldType::view_t, typename DetectorType::algebra_type,
                         detray::constrained_step<>>;
  using NavigatorType = detray::navigator<const DetectorType>;
  using FitterType = traccc::kalman_fitter<StepperType, NavigatorType>;

  using SiliconDetectorDescriptionType = traccc::silicon_detector_description;

  using ClusterizationAlgorithmType = traccc::host::clusterization_algorithm;
  using SpacepointFormationAlgorithmType =
      traccc::host::spacepoint_formation_algorithm;
  using SeedingAlgorithmType = traccc::seeding_algorithm;
  using TrackParametersEstimationAlgorithmType =
      traccc::track_params_estimation;
  using FindingAlgorithmType =
      traccc::finding_algorithm<StepperType, NavigatorType>;
  using FittingAlgorithmType = traccc::fitting_algorithm<FitterType>;
  using AmbiguityResolutionAlgorithmType =
      traccc::greedy_ambiguity_resolution_algorithm;
};

}  // namespace ActsExamples::Traccc::Host
