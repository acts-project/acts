#include <iostream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/GSLStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Observers.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"

using namespace Acts;

int
main()
{
  typedef ConstantFieldSvc         BField_type;
  typedef GSLStepper<BField_type>  Stepper_type;
  typedef Propagator<Stepper_type> Propagator_type;

  BField_type::Config c;
  c.field = {0, 0, 2 * units::_T};
  BField_type     bField(std::move(c));
  Stepper_type    gsl_stepper(std::move(bField));
  Propagator_type propagator(std::move(gsl_stepper));

  typedef ObserverList<PathLengthObserver> ObsList_type;
  typedef AbortList<MaxPathLength>         AbortConditions_type;

  Propagator_type::Options<ObsList_type, AbortConditions_type> options;
  AbortConditions_type& al              = options.stop_conditions;
  al.get<MaxPathLength>().maxPathLength = 5 * units::_m;

  Vector3D              pos(0, 0, 0);
  Vector3D              mom(1 * units::_GeV, 0, 0);
  CurvilinearParameters pars(nullptr, pos, mom, +1);
  double                totalPathLength = 0;
  for (unsigned int i = 0; i < 10000; ++i) {
    auto r = propagator.propagate(pars, options);
    totalPathLength += r.get<PathLengthObserver::result_type>().pathLength;
  }

  std::cout << totalPathLength / units::_cm << std::endl;
  return 0;
}
