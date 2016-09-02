#include "ACTS/Extrapolation/Propagator.hpp"
#include <random>
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
  // random numbers for smearing measurements
  std::default_random_engine             e;
  std::uniform_real_distribution<double> uniform(-1, 1);

  typedef ConstantFieldSvc      BField_t;
  typedef GSLStepper<BField_t>  Stepper_t;
  typedef Propagator<Stepper_t> Propagator_t;

  BField_t::Config c;
  c.field = {0, 0, 2 * units::_T};
  BField_t     bField(std::move(c));
  Stepper_t    gsl_stepper(std::move(bField));
  Propagator_t propagator(std::move(gsl_stepper));

  typedef Propagator_t::observer_list_t<CurvilinearParameters,
                                        PathLengthObserver,
                                        HitSimulator>
      ObsList_t;

  ObsList_t     ol;
  HitSimulator& hit_sim = ol.get<HitSimulator>();
  typedef std::tuple<float, float, float> tuple_t;
  hit_sim.barrel
      = {tuple_t(5 * units::_cm, -50 * units::_cm, 50 * units::_cm),
         tuple_t(15 * units::_cm, -55 * units::_cm, 55 * units::_cm),
         tuple_t(25 * units::_cm, -60 * units::_cm, 60 * units::_cm),
         tuple_t(50 * units::_cm, -75 * units::_cm, 75 * units::_cm),
         tuple_t(60 * units::_cm, -80 * units::_cm, 80 * units::_cm)};

  hit_sim.endcaps
      = {tuple_t(70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
         tuple_t(85 * units::_cm, 15 * units::_cm, 60 * units::_cm),
         tuple_t(-70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
         tuple_t(-85 * units::_cm, 15 * units::_cm, 60 * units::_cm)};

  HitSimulator::result_type hits;
  double                    totalPathLength = 0;
  for (unsigned int i = 0; i < 1000; ++i) {
    Vector3D pos(0, 0, 0);
    Vector3D mom(uniform(e), uniform(e), uniform(e));
    mom *= 1 * units::_GeV;
    CurvilinearParameters pars(nullptr, pos, mom, +1);
    auto                  r = propagator.propagate(pars, ol);
    totalPathLength += r.get<PathLengthObserver::result_type>().pathLength;
    hits.add(r.get<HitSimulator::result_type>());
  }
  std::cout << "total path length = " << totalPathLength / units::_m
            << std::endl;
  std::cout << "simulated hits:" << std::endl;
  hits.print();

  return 0;
}
