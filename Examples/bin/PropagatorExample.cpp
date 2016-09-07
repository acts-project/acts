#include "ACTS/Extrapolation/Propagator.hpp"
#include <cstdio>
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
  typedef ConstantFieldSvc         BField_type;
  typedef GSLStepper<BField_type>  Stepper_type;
  typedef Propagator<Stepper_type> Propagator_type;

  BField_type::Config c;
  c.field = {0, 0, 2 * units::_T};
  BField_type     bField(std::move(c));
  Stepper_type    gsl_stepper(std::move(bField));
  Propagator_type propagator(std::move(gsl_stepper));

  // first example: simple propagation of a single track
  {
    typedef ObserverList<PathLengthObserver, DebugObserver> ObsList_type;
    typedef AbortList<MaxPathLength, MaxRadius> AbortConditions_type;

    Propagator_type::Options<ObsList_type, AbortConditions_type> options;
    AbortConditions_type& al              = options.stop_conditions;
    al.get<MaxPathLength>().maxPathLength = 35 * units::_cm;
    al.get<MaxRadius>().maxRadius         = 15 * units::_cm;

    std::FILE*    file          = std::fopen("track.txt", "w");
    ObsList_type& ol            = options.observer_list;
    ol.get<DebugObserver>().out = file;
    Vector3D              pos(0, 0, 0);
    Vector3D              mom(1 * units::_GeV, 0, 0);
    CurvilinearParameters pars(nullptr, pos, mom, +1);
    auto                  r = propagator.propagate(pars, options);
    std::cout << "path length = "
              << r.get<PathLengthObserver::result_type>().pathLength / units::_m
              << std::endl;
    std::fclose(file);
  }

  // second example: use propagator to generate pseudo-hits
  {
    // random numbers for initializing track parameters
    std::random_device                     r;
    std::default_random_engine             e(r());
    std::uniform_real_distribution<double> uniform(0, 1);

    typedef ObserverList<HitSimulator> ObsList_type;
    typedef AbortList<>                AbortConditions_type;

    Propagator_type::Options<ObsList_type, AbortConditions_type> options;
    ObsList_type& ol      = options.observer_list;
    HitSimulator& hit_sim = ol.get<HitSimulator>();
    typedef std::tuple<float, float, float> tuple_type;

    // setup barrel structure: r, zmin, zmax
    hit_sim.barrel
        = {tuple_type(5 * units::_cm, -50 * units::_cm, 50 * units::_cm),
           tuple_type(15 * units::_cm, -55 * units::_cm, 55 * units::_cm),
           tuple_type(25 * units::_cm, -60 * units::_cm, 60 * units::_cm),
           tuple_type(50 * units::_cm, -75 * units::_cm, 75 * units::_cm),
           tuple_type(60 * units::_cm, -80 * units::_cm, 80 * units::_cm)};

    // setup endcap structure: z, rmin, rmax
    hit_sim.endcaps
        = {tuple_type(70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
           tuple_type(85 * units::_cm, 15 * units::_cm, 60 * units::_cm),
           tuple_type(-70 * units::_cm, 5 * units::_cm, 45 * units::_cm),
           tuple_type(-85 * units::_cm, 15 * units::_cm, 60 * units::_cm)};

    HitSimulator::result_type hits;
    for (unsigned int i = 0; i < 1000; ++i) {
      Vector3D pos(0, 0, 0);
      double   pT  = 10 * units::_GeV;
      double   eta = -4 + uniform(e) * 8;
      double   phi = -M_PI + uniform(e) * 2 * M_PI;
      Vector3D mom(pT * cos(phi), pT * sin(phi), pT / tan(2 * atan(exp(-eta))));
      mom *= 1 * units::_GeV;
      CurvilinearParameters pars(nullptr, pos, mom, +1);
      auto                  r = propagator.propagate(pars, options);
      hits.add(r.get<HitSimulator::result_type>());
    }

    std::ofstream hits_file("hits.txt");
    hits.print(hits_file);
  }

  return 0;
}
