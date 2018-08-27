// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/program_options.hpp>
#include <iostream>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Extrapolation/RungeKuttaEngine.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace po = boost::program_options;
using namespace Acts;

int
main(int argc, char* argv[])
{
  unsigned int toys    = 1;
  double       pT      = 1;
  double       Bz      = 1;
  double       maxPath = 1;
  unsigned int lvl     = Acts::Logging::INFO;
  bool         withCov = true;

  try {
    po::options_description desc("Allowed options");
    // clang-format off
  desc.add_options()
      ("help", "produce help message")
      ("toys",po::value<unsigned int>(&toys)->default_value(10000),"number of tracks to propagate")
      ("pT",po::value<double>(&pT)->default_value(1),"transverse momentum in GeV")
      ("B",po::value<double>(&Bz)->default_value(2),"z-component of B-field in T")
      ("path",po::value<double>(&maxPath)->default_value(5),"maximum path length in m")
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

  auto myLogger = getDefaultLogger("ATLAS_Stepper", Acts::Logging::Level(lvl));
  ACTS_LOCAL_LOGGER(std::move(myLogger));

  // print information about profiling setup
  ACTS_INFO("propagating " << toys << " tracks with pT = " << pT << "GeV in a "
                           << Bz
                           << "T B-field");

  using BField_type   = ConstantBField;
  auto magnetic_field = std::make_unique<const BField_type>(0, 0, Bz / 1000.);

  RungeKuttaEngine<>::Config c;
  c.fieldService = std::move(magnetic_field);
  c.pathLimit    = maxPath * units::_m;

  RungeKuttaEngine<> propagator(c);

  CylinderSurface surface(nullptr, 100 * units::_m, 30 * units::_m);

  Vector3D          pos(0, 0, 0);
  Vector3D          mom(pT * units::_GeV, 0, 0);
  ActsSymMatrixD<5> cov;
  cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);

  std::unique_ptr<const ActsSymMatrixD<5>> covPtr = nullptr;
  if (withCov) {
    covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
  }
  CurvilinearParameters pars(std::move(covPtr), pos, mom, +1);

  ExtrapolationCell<TrackParameters> exCell(pars);
  exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);

  double totalPathLength = 0;
  for (unsigned int i = 0; i < toys; ++i) {
    exCell.pathLength         = 0;
    exCell.leadParameters     = &pars;
    exCell.lastLeadParameters = 0;
    propagator.propagate(exCell, surface);
    ACTS_DEBUG("reached position (" << exCell.leadParameters->position().x()
                                    << ", "
                                    << exCell.leadParameters->position().y()
                                    << ", "
                                    << exCell.leadParameters->position().z()
                                    << ") in "
                                    << exCell.nSteps
                                    << " steps");
    totalPathLength += exCell.pathLength;
  }

  ACTS_INFO("average path length = " << totalPathLength / toys / units::_mm
                                     << "mm");

  return 0;
}
