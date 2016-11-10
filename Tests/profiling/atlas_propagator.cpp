#include <boost/program_options.hpp>
#include <iostream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Logger.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace po = boost::program_options;
using namespace Acts;

int
main(int argc, char* argv[])
{
  unsigned int toys    = 1;
  double       pT      = 1;
  double       Bz      = 1;
  double       maxPath = 1;

  try {
    po::options_description desc("Allowed options");
    // clang-format off
  desc.add_options()
      ("help", "produce help message")
      ("toys",po::value<unsigned int>(&toys)->default_value(10000),"number of tracks to propagate")
      ("pT",po::value<double>(&pT)->default_value(1 * units::_GeV),"transverse momentum in GeV")
      ("B",po::value<double>(&Bz)->default_value(2),"z-component of B-field in T")
      ("path",po::value<double>(&maxPath)->default_value(5 * units::_m),"maximum path length in m");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 0;
    }
  } catch (std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }

  // print information about profiling setup
  std::cout << "propagating " << toys
            << " tracks with pT = " << pT / units::_GeV << "GeV in a " << Bz
            << "T B-field" << std::endl;

  typedef ConstantBField       BField_type;
  std::unique_ptr<BField_type> magnetic_field
      = std::make_unique<BField_type>(0, 0, Bz);

  RungeKuttaEngine<>::Config c;
  c.fieldService  = std::move(magnetic_field);
  c.maxPathLength = maxPath;

  RungeKuttaEngine<> propagator(c);

  CylinderSurface surface(nullptr, 100 * units::_m, 30 * units::_m);

  Vector3D          pos(0, 0, 0);
  Vector3D          mom(pT, 0, 0);
  ActsSymMatrixD<5> cov;
  cov << 10 * units::_mm, 0, 0, 0, 0, 0, 10 * units::_mm, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1. / (10 * units::_GeV);

  CurvilinearParameters pars(
      std::make_unique<ActsSymMatrixD<5>>(cov), pos, mom, +1);
  ExtrapolationCell<TrackParameters> exCell(pars);
  exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);

  double totalPathLength = 0;
  for (unsigned int i = 0; i < toys; ++i) {
    exCell.pathLength         = 0;
    exCell.leadParameters     = &pars;
    exCell.lastLeadParameters = 0;
    propagator.propagate(exCell, surface);
    totalPathLength += exCell.pathLength;
  }

  std::cout << "average path length = " << totalPathLength / toys / units::_mm
            << "mm" << std::endl;
  return 0;
}
