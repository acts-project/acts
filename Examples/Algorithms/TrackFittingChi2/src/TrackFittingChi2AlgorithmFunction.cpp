// This file is part of the Acts project.
//
// Copyright (C) 2021-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Direction.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/TrackFitting/Chi2Fitter.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/TrackFittingChi2/TrackFittingChi2Algorithm.hpp"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
class SourceLink;
class TrackingGeometry;
class VectorMultiTrajectory;
}  // namespace Acts

namespace {

using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter =
    Acts::Experimental::Chi2Fitter<Propagator, Acts::VectorMultiTrajectory>;

struct TrackFitterChi2FunctionImpl
    : public ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Function {
  Fitter trackFitterChi2;

  TrackFitterChi2FunctionImpl(Fitter&& f) : trackFitterChi2(std::move(f)) {}

  ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Result operator()(
      const std::vector<Acts::SourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Options&
          options,
      ActsExamples::TrackContainer& trackContainer) const override {
    return trackFitterChi2.fit(sourceLinks.begin(), sourceLinks.end(),
                               initialParameters, options, trackContainer);
  };
};

}  // namespace

std::shared_ptr<
    ActsExamples::TrackFittingChi2Algorithm::TrackFitterChi2Function>
ActsExamples::TrackFittingChi2Algorithm::makeTrackFitterChi2Function(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField) {
  Stepper stepper(std::move(magneticField));
  Acts::Navigator::Config cfg{std::move(trackingGeometry)};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitterChi2(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  return std::make_shared<TrackFitterChi2FunctionImpl>(
      std::move(trackFitterChi2));
}
