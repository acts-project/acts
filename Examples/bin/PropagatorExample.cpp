#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/GSLStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Observers.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"

using namespace Acts;

int
main()
{
  typedef ConstantFieldSvc      BField_t;
  typedef GSLStepper<BField_t>  Stepper_t;
  typedef Propagator<Stepper_t> Propagator_t;

  BField_t::Config c;
  c.field = {0, 0, 2 * units::_T};
  BField_t     bField(std::move(c));
  Stepper_t    gsl_stepper(std::move(bField));
  Propagator_t propagator(std::move(gsl_stepper));

  Vector3D              pos(0, 0, 0);
  Vector3D              mom(1 * units::_GeV, 0, 0);
  CurvilinearParameters pars(nullptr, pos, mom, +1);

  typedef Propagator_t::observer_list_t<CurvilinearParameters,
                                        DebugObserver,
                                        PathLengthObserver>
       ObsList_t;
  auto r = propagator.propagate(pars, ObsList_t());

  std::cout << "total path length = "
            << r.get<PathLengthObserver::result_type>().pathLength / units::_m
            << std::endl;

  return 0;
}
