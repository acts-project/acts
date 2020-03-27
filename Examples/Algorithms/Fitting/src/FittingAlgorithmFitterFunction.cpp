// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <iostream>
#include <map>
#include <random>
#include <stdexcept>

#include "ACTFW/Fitting/FittingAlgorithm.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace {
template <typename Fitter>
struct FitterFunctionImpl {
  Fitter fitter;

  FitterFunctionImpl(Fitter&& f) : fitter(std::move(f)) {}

  FW::FittingAlgorithm::FitterResult operator()(
      const std::vector<FW::SimSourceLink>& sourceLinks,
      const FW::TrackParameters& initialParameters,
      const Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>& options) const {
    return fitter.fit(sourceLinks, initialParameters, options);
  };
};
}  // namespace

FW::FittingAlgorithm::FitterFunction FW::FittingAlgorithm::makeFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    Options::BFieldVariant magneticField, Acts::Logging::Level lvl) {
  using Updater = Acts::GainMatrixUpdater<Acts::BoundParameters>;
  using Smoother = Acts::GainMatrixSmoother<Acts::BoundParameters>;

  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [trackingGeometry, lvl](auto&& inputField) -> FitterFunction {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        using Stepper = Acts::EigenStepper<MagneticField>;
        using Navigator = Acts::Navigator;
        using Propagator = Acts::Propagator<Stepper, Navigator>;
        using Fitter = Acts::KalmanFitter<Propagator, Updater, Smoother>;

        // construct all components for the fitter
        MagneticField field(std::move(inputField));
        Stepper stepper(std::move(field));
        Navigator navigator(trackingGeometry);
        navigator.resolvePassive = false;
        navigator.resolveMaterial = true;
        navigator.resolveSensitive = true;
        Propagator propagator(std::move(stepper), std::move(navigator));
        Fitter fitter(std::move(propagator),
                      Acts::getDefaultLogger("KalmanFitter", lvl));

        // build the fitter functions. owns the fitter object.
        return FitterFunctionImpl<Fitter>(std::move(fitter));
      },
      std::move(magneticField));
}
