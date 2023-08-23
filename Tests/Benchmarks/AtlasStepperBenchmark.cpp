// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/AtlasStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <iostream>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace Acts;
using namespace Acts::UnitLiterals;

int main(int argc, char* argv[]) {
  unsigned int toys = 1;
  double ptInGeV = 1;
  double BzInT = 1;
  double maxPathInM = 1;
  unsigned int lvl = Acts::Logging::INFO;
  bool withCov = true;

  // Create a test context
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();

  try {
    po::options_description desc("Allowed options");
    // clang-format off
  desc.add_options()
      ("help", "produce help message")
      ("toys",po::value<unsigned int>(&toys)->default_value(20000),"number of tracks to propagate")
      ("pT",po::value<double>(&ptInGeV)->default_value(1),"transverse momentum in GeV")
      ("B",po::value<double>(&BzInT)->default_value(2),"z-component of B-field in T")
      ("path",po::value<double>(&maxPathInM)->default_value(5),"maximum path length in m")
      ("cov",po::value<bool>(&withCov)->default_value(true),"propagation with covariance matrix")
      ("verbose",po::value<unsigned int>(&lvl)->default_value(Acts::Logging::INFO),"logging level");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help") != 0u) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  ACTS_LOCAL_LOGGER(
      getDefaultLogger("ATLAS_Stepper", Acts::Logging::Level(lvl)));

  // print information about profiling setup
  ACTS_INFO("propagating " << toys << " tracks with pT = " << ptInGeV
                           << "GeV in a " << BzInT << "T B-field");

  using BField_type = ConstantBField;
  using Stepper_type = AtlasStepper;
  using Propagator_type = Propagator<Stepper_type>;
  using Covariance = BoundSquareMatrix;

  auto bField =
      std::make_shared<BField_type>(Vector3{0, 0, BzInT * UnitConstants::T});
  Stepper_type atlas_stepper(std::move(bField));
  Propagator_type propagator(std::move(atlas_stepper));

  PropagatorOptions<> options(tgContext, mfContext);
  options.pathLimit = maxPathInM * UnitConstants::m;

  Covariance cov;
  // clang-format off
  cov << 10_mm, 0, 0, 0, 0, 0,
         0, 10_mm, 0, 0, 0, 0,
         0, 0, 1, 0, 0, 0,
         0, 0, 0, 1, 0, 0,
         0, 0, 0, 0, 1_e / 10_GeV, 0,
         0, 0, 0, 0, 0, 0;
  // clang-format on

  std::optional<Covariance> optCov = std::nullopt;
  if (withCov) {
    optCov = cov;
  }
  CurvilinearTrackParameters pars(Vector4::Zero(), 0_degree, 90_degree,
                                  1_e / ptInGeV * UnitConstants::GeV, optCov,
                                  ParticleHypothesis::pion());

  double totalPathLength = 0;
  size_t num_iters = 0;
  const auto propagation_bench_result = Acts::Test::microBenchmark(
      [&] {
        auto r = propagator.propagate(pars, options).value();
        if (totalPathLength == 0.) {
          ACTS_DEBUG("reached position "
                     << r.endParameters->position(tgContext).transpose()
                     << " in " << r.steps << " steps");
        }
        totalPathLength += r.pathLength;
        ++num_iters;
        return r;
      },
      1, toys);

  ACTS_INFO("Execution stats: " << propagation_bench_result);
  ACTS_INFO("average path length = " << totalPathLength / num_iters / 1_mm
                                     << "mm");
  return 0;
}
