#include <boost/program_options.hpp>
#include <iostream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/AtlasStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Observers.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace po = boost::program_options;
using namespace Acts;

int
main(int argc, char* argv[])
{
  int    toys    = 1;
  double pT      = 1;
  double Bz      = 1;
  double maxPath = 1;

  try {
    po::options_description desc("Allowed options");
    // clang-format off
  desc.add_options()
      ("help", "produce help message")
      ("toys",po::value<int>(&toys)->default_value(10000),"number of tracks to propagate")
      ("pT",po::value<double>(&pT)->default_value(1 * units::_GeV),"transverse momentum in GeV")
      ("B",po::value<double>(&Bz)->default_value(2 * units::_T),"z-component of B-field in T")
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
            << " tracks with pT = " << pT / units::_GeV << "GeV in a "
            << Bz / units::_T << "T B-field" << std::endl;

  typedef ConstantFieldSvc          BField_type;
  typedef AtlasStepper<BField_type> Stepper_type;
  typedef Propagator<Stepper_type>  Propagator_type;

  BField_type::Config c;
  c.field = {0, 0, Bz};
  BField_type     bField(std::move(c));
  Stepper_type    atlas_stepper(std::move(bField));
  Propagator_type propagator(std::move(atlas_stepper));

  typedef ObserverList<PathLengthObserver> ObsList_type;
  typedef AbortList<MaxPathLength>         AbortConditions_type;

  Propagator_type::Options<forward, ObsList_type, AbortConditions_type> options;
  AbortConditions_type& al              = options.stop_conditions;
  al.get<MaxPathLength>().maxPathLength = maxPath;

  Vector3D              pos(0, 0, 0);
  Vector3D              mom(pT, 0, 0);
  CurvilinearParameters pars(nullptr, pos, mom, +1);
  double                totalPathLength = 0;
  for (unsigned int i = 0; i < 10000; ++i) {
    auto r = propagator.propagate(pars, options);
    totalPathLength += r.get<PathLengthObserver::result_type>().pathLength;
  }

  std::cout << "average path length = " << totalPathLength / toys / units::_mm
            << "mm" << std::endl;
  return 0;
}
