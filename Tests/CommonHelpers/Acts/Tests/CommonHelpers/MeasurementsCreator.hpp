// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <random>
#include <vector>

namespace Acts {
namespace Test {

/// All supported simulated measurement types.
enum class MeasurementType {
  eLoc0,
  eLoc1,
  eLoc01,
};

/// Measurement resolution configuration for a single detector type.
struct MeasurementResolution {
  MeasurementType type = MeasurementType::eLoc0;
  // Depending on the type, only the first value is used.
  std::array<double, 2> stddev = {50 * UnitConstants::um, 0};
};

/// Measurement resolution configuration for a full detector geometry.
using MeasurementResolutionMap =
    Acts::GeometryHierarchyMap<MeasurementResolution>;

/// Result struct for generated measurements and outliers.
struct Measurements {
  // use shared_ptr to provide stable memory addresses. unique_ptr would be
  // clearer, but that makes the type non-copyable and very hard to use.
  std::vector<std::shared_ptr<Acts::FittableMeasurement<TestSourceLink>>> store;
  std::vector<TestSourceLink> sourceLinks;
  std::vector<TestSourceLink> outlierSourceLinks;
};

/// Propagator action to create smeared measurements.
struct MeasurementsCreator {
  using result_type = Measurements;

  MeasurementResolutionMap resolutions;
  std::default_random_engine* rng = nullptr;
  size_t sourceId = 0;
  // how far away from the measurements the outliers should be
  double distanceOutlier = 10 * Acts::UnitConstants::mm;

  /// @brief Operater that is callable by an ActionList. The function
  /// collects the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] result Vector of matching surfaces
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    using namespace Acts::UnitLiterals;
    using VariantMeasurement = Acts::FittableMeasurement<TestSourceLink>;
    using MeasurementLoc0 =
        Acts::Measurement<TestSourceLink, Acts::BoundIndices, Acts::eBoundLoc0>;
    using MeasurementLoc1 =
        Acts::Measurement<TestSourceLink, Acts::BoundIndices, Acts::eBoundLoc1>;
    using MeasurementLoc01 =
        Acts::Measurement<TestSourceLink, Acts::BoundIndices, Acts::eBoundLoc0,
                          Acts::eBoundLoc1>;

    const auto& logger = state.options.logger;

    // only generate measurements on surfaces
    if (not state.navigation.currentSurface) {
      return;
    }
    const Acts::Surface& surface = *state.navigation.currentSurface;
    const Acts::GeometryIdentifier geoId = surface.geometryId();
    // only generate measurements on sensitive surface
    if (not geoId.sensitive()) {
      ACTS_VERBOSE("Create no measurements on non-sensitive surface " << geoId);
      return;
    }
    // only generate measurements if a resolution is configured
    auto found = resolutions.find(geoId);
    if (found == resolutions.end()) {
      ACTS_VERBOSE("No resolution configured for sensitive surface " << geoId);
      return;
    }
    const MeasurementResolution& resolution = *found;

    // Apply global to local
    Acts::Vector2D loc =
        surface
            .globalToLocal(state.geoContext, stepper.position(state.stepping),
                           stepper.direction(state.stepping))
            .value();

    std::normal_distribution<double> normalDist(0., 1.);

    // compute covariance for all components, might contain bogus values
    // depending on the configuration. but those remain unused.
    Vector2D stddev(resolution.stddev[0], resolution.stddev[1]);
    SymMatrix2D cov = stddev.cwiseProduct(stddev).asDiagonal();

    const VariantMeasurement* measPtr = nullptr;
    const VariantMeasurement* outlierPtr = nullptr;
    if (resolution.type == MeasurementType::eLoc0) {
      // create measurement w/ dummy source link
      double parMeas = loc[0] + stddev[0] * normalDist(*rng);
      measPtr = result.store
                    .emplace_back(std::make_shared<VariantMeasurement>(
                        MeasurementLoc0(surface.getSharedPtr(), {},
                                        cov.topLeftCorner<1, 1>(), parMeas)))
                    .get();
      // create outlier w/ dummy source link
      double parOutlier = parMeas + distanceOutlier;
      outlierPtr =
          result.store
              .emplace_back(std::make_shared<VariantMeasurement>(
                  MeasurementLoc0(surface.getSharedPtr(), {},
                                  cov.topLeftCorner<1, 1>(), parOutlier)))
              .get();
    } else if (resolution.type == MeasurementType::eLoc1) {
      // create measurement w/ dummy source link
      // yes, using stddev[0] and cov(0,0) is correct here. this accesses the
      // first configuration parameter not the first local coordinate.
      double parMeas = loc[1] + stddev[0] * normalDist(*rng);
      measPtr = result.store
                    .emplace_back(std::make_shared<VariantMeasurement>(
                        MeasurementLoc1(surface.getSharedPtr(), {},
                                        cov.topLeftCorner<1, 1>(), parMeas)))
                    .get();
      // create outlier w/ dummy source link
      double parOutlier = parMeas + distanceOutlier;
      outlierPtr =
          result.store
              .emplace_back(std::make_shared<VariantMeasurement>(
                  MeasurementLoc1(surface.getSharedPtr(), {},
                                  cov.topLeftCorner<1, 1>(), parOutlier)))
              .get();
    } else if (resolution.type == MeasurementType::eLoc01) {
      // create measurement w/ dummy source link
      Vector2D parMeas = loc;
      loc[0] += stddev[0] * normalDist(*rng);
      loc[1] += stddev[1] * normalDist(*rng);
      measPtr =
          result.store
              .emplace_back(std::make_shared<VariantMeasurement>(
                  MeasurementLoc01(surface.getSharedPtr(), {}, cov, parMeas)))
              .get();
      // create outlier w/ dummy source link
      Vector2D parOutlier = parMeas;
      parOutlier[0] += distanceOutlier;
      parOutlier[1] -= distanceOutlier;
      outlierPtr = result.store
                       .emplace_back(std::make_shared<VariantMeasurement>(
                           MeasurementLoc01(surface.getSharedPtr(), {}, cov,
                                            parOutlier)))
                       .get();
    }

    result.sourceLinks.emplace_back(*measPtr, sourceId);
    result.outlierSourceLinks.emplace_back(*outlierPtr, sourceId);
  }
};

/// Propagate the track create smeared measurements from local coordinates.
template <typename propagator_t, typename track_parameters_t>
Measurements createMeasurements(const propagator_t& propagator,
                                const Acts::GeometryContext& geoCtx,
                                const Acts::MagneticFieldContext& magCtx,
                                const track_parameters_t& trackParameters,
                                const MeasurementResolutionMap& resolutions,
                                std::default_random_engine& rng,
                                size_t sourceId = 0u) {
  using Actions = Acts::ActionList<MeasurementsCreator>;
  using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;

  // Set options for propagator
  auto logger =
      Acts::getDefaultLogger("MeasurementCreator", Acts::Logging::INFO);
  Acts::PropagatorOptions<Actions, Aborters> options(
      geoCtx, magCtx, Acts::LoggerWrapper(*logger));
  auto& creator = options.actionList.get<MeasurementsCreator>();
  creator.resolutions = resolutions;
  creator.rng = &rng;
  creator.sourceId = sourceId;

  // Launch and collect the measurements
  return propagator.propagate(trackParameters, options)
      .value()
      .template get<Measurements>();
}

}  // namespace Test
}  // namespace Acts
