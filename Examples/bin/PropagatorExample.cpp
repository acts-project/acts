#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/GSLStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"

using namespace Acts;

int
main()
{
  typedef FWE::ConstantFieldSvc BField_t;
  typedef GSLStepper<BField_t>  Stepper_t;
  BField_t::Config              c;
  c.field = {0, 0, 2};
  BField_t              bField(std::move(c));
  Stepper_t             gsl_stepper(std::move(bField));
  Propagator<Stepper_t> propagator(std::move(gsl_stepper));

  Vector3D              pos(0, 0, 0);
  Vector3D              mom(1, 0, 0);
  CurvilinearParameters pars(nullptr, pos, mom, +1);

  std::cout << pars << std::endl;
  auto r = propagator.propagate(pars);
  std::cout << r.endParameters << std::endl;

  return 0;
}
