// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "ActsExamples/Fitting/DirectedFittingAlgorithm.hpp"
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"

#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

#include <boost/program_options.hpp>

namespace {
template <typename Fitter>
struct FitterFunctionImpl {
  Fitter fitter;

  FitterFunctionImpl(Fitter&& f) : fitter(std::move(f)) {}

  ActsExamples::DirectedFittingAlgorithm::FitterResult operator()(
      const std::vector<ActsExamples::SimSourceLink>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& options,
      const std::vector<const Acts::Surface*>& sSequence) const {
    return fitter.fit(sourceLinks, initialParameters, options, sSequence);
  };
};
}  // namespace

ActsExamples::DirectedFittingAlgorithm::FitterFunction
ActsExamples::DirectedFittingAlgorithm::makeFitterFunction(
    Options::BFieldVariant magneticField) {
  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [](auto&& inputField) -> FitterFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<MagneticField>;
        using Navigator = Acts::DirectNavigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

        // construct all components for the fitter
        MagneticField field(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter fitter(std::move(propagator));

        // build the fitter functions. owns the fitter object.
        return FitterFunctionImpl<Fitter>(std::move(fitter));
      },
      std::move(magneticField));
}
