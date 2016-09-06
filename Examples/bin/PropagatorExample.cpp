#include "ACTS/Extrapolation/Propagator.hpp"
#include <fstream>
#include <random>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
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

  // first example: simple propagation of a single track
  {
    typedef ObserverList<PathLengthObserver, DebugObserver> ObsList_t;
    typedef AbortList<MaxPathLength, MaxRadius>             AbortConditions_t;
    AbortConditions_t al;
    al.get<MaxPathLength>().maxPathLength = 35 * units::_cm;
    al.get<MaxRadius>().maxRadius         = 15 * units::_cm;

    std::ofstream out_file("track.txt");
    ObsList_t     ol;
    ol.get<DebugObserver>().out.rdbuf(out_file.rdbuf());
    Vector3D              pos(0, 0, 0);
    Vector3D              mom(1 * units::_GeV, 0, 0);
    CurvilinearParameters pars(nullptr, pos, mom, +1);
    auto                  r = propagator.propagate(pars, ol, al);
    std::cout << "path length = "
              << r.get<PathLengthObserver::result_type>().pathLength / units::_m
              << std::endl;
  }

  // second example: use propagator to generate pseudo-hits
  {
    // random numbers for initializing track parameters
    std::random_device                     r;
    std::default_random_engine             e(r());
    std::uniform_real_distribution<double> uniform(0, 1);

    typedef ObserverList<HitSimulator> ObsList_t;
    typedef AbortList<>                AbortConditions_t;

    ObsList_t     ol;
    HitSimulator& hit_sim = ol.get<HitSimulator>();
    typedef std::tuple<float, float, float> tuple_t;

    // setup barrel structure: r, zmin, zmax
    hit_sim.barrel
        = {tuple_t(5 * units::_cm, -50 * units::_cm, 50 * units::_cm),
           tuple_t(15 * units::_cm, -55 * units::_cm, 55 * units::_cm),
           tuple_t(25 * units::_cm, -60 * units::_cm, 60 * units::_cm),
           tuple_t(50 * units::_cm, -75 * units::_cm, 75 * units::_cm),
           tuple_t(60 * units::_cm, -80 * units::_cm, 80 * units::_cm)};

    // setup endcap structure: z, rmin, rmax
    hit_sim.endcaps
        = {tuple_t(70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
           tuple_t(85 * units::_cm, 15 * units::_cm, 60 * units::_cm),
           tuple_t(-70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
           tuple_t(-85 * units::_cm, 15 * units::_cm, 60 * units::_cm)};

    HitSimulator::result_type hits;
    for (unsigned int i = 0; i < 1000; ++i) {
      Vector3D pos(0, 0, 0);
      double   pT  = 10 * units::_GeV;
      double   eta = -4 + uniform(e) * 8;
      double   phi = -M_PI + uniform(e) * 2 * M_PI;
      Vector3D mom(pT * cos(phi), pT * sin(phi), pT / tan(2 * atan(exp(-eta))));
      mom *= 1 * units::_GeV;
      CurvilinearParameters pars(nullptr, pos, mom, +1);
      auto r = propagator.propagate(pars, ol, AbortConditions_t());
      hits.add(r.get<HitSimulator::result_type>());
    }

    std::ofstream hits_file("hits.txt");
    hits.print(hits_file);
  }

  return 0;
}
