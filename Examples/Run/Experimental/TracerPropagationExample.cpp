// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Experimental/Tracer.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <random>

#include "BuildTracerDetector.hpp"

Acts::BoundTrackParameters generateParameters(
    std::shared_ptr<const Acts::PerigeeSurface> psf,
    std::uniform_real_distribution<double>& phiDist,
    std::uniform_real_distribution<double>& etaDist, std::mt19937& gen) {
  double phi = phiDist(gen);
  double eta = etaDist(gen);
  double theta = 2 * std::atan(std::exp(-eta));

  Acts::BoundVector pars;
  pars << 0., 0., phi, theta, 0.001, 0.;

  // charged extrapolation - with hit recording
  return Acts::BoundTrackParameters(psf, std::move(pars), std::nullopt);
}

int main(int argc, char* argv[]) {
  const std::string jsonFile = argv[1];

  const std::string outputFile = argv[2];
  bool writeOut = outputFile != "none";

  auto detector = Acts::detectorFromJson(jsonFile);

  Acts::Tracer::Config tracerCfg;
  tracerCfg.world = detector;
  Acts::Tracer tracer(tracerCfg);

  Acts::StraightLineStepper slStepper;

  Acts::Propagator<Acts::StraightLineStepper, Acts::Tracer> propagator(
      std::move(slStepper), std::move(tracer));

  const int nEvents = std::atoi(argv[3]);
  const double etaValue = std::atof(argv[4]);

  bool targetEta = argc > 5;

  double minEta = targetEta ? (double)etaValue : -std::abs(double(etaValue));
  double maxEta = double(etaValue);

  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());

  // Setup random number distributions for some quantities
  std::uniform_real_distribution<double> phiDist(-M_PI, M_PI);
  std::uniform_real_distribution<double> etaDist(minEta, maxEta);

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3(0., 0., 0.));

  Acts::GeometryContext gctx{};
  Acts::MagneticFieldContext mctx{};

  auto exampleLogger =  Acts::getDefaultLogger("TracerExample", Acts::Logging::INFO);

  Acts::ActionList<Acts::detail::SteppingLogger> actionList;
  auto& steppingLogger = actionList.get<Acts::detail::SteppingLogger>();
  steppingLogger.sterile = not writeOut;

  Acts::PropagatorOptions<decltype(actionList)> pOptions(gctx, mctx, Acts::LoggerWrapper{*exampleLogger});

  std::ofstream output;
  if (writeOut){
     output.open(outputFile);
  }
 
  for (int ie = 0; ie < nEvents; ++ie) {
    auto parameters = generateParameters(surface, phiDist, etaDist, gen);
    auto propagateResult = propagator.propagate(parameters, pOptions);
    if (propagateResult.ok() and writeOut){
       auto steppingResults =
            propagateResult.value().get<Acts::detail::SteppingLogger::result_type>();
       for (const auto& s : steppingResults.steps){
         std::string positionStr = Acts::toString(s.position);
         positionStr.erase(0u,1);
         positionStr.erase(positionStr.size()-1,1);
         output << ie << ", " << positionStr << '\n';
       }
    }
  }

  if (writeOut){
    output.close();
  }
}
