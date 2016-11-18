#include "ACTS/Extrapolation/Propagator.hpp"
#include <cstdio>
#include <fstream>
#include <list>
#include <random>
#include <tuple>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/AtlasStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Utilities/Units.hpp"

using namespace Acts;
using namespace Acts::propagation;

struct DebugObserver
{
  template <typename TrackParameters>
  void
  operator()(const TrackParameters& current,
             const TrackParameters& /*previous*/) const
  {
    const auto& pos = current.position();
    const auto& mom = current.momentum();
    fprintf(out,
            "x = %.2fmm, y = %.2fmm, z = %.2fmm, pT = %.3fGeV, phi = %.5f, "
            "theta = %.5f\n",
            pos(0) / units::_mm,
            pos(1) / units::_mm,
            pos(2) / units::_mm,
            mom.perp() / units::_GeV,
            mom.phi(),
            mom.theta());
  }

  std::FILE* out = stdout;
};

struct HitSimulator
{
  struct this_result
  {
    std::list<Vector3D> hits = {};

    void
    add(const this_result& other)
    {
      hits.insert(hits.end(), other.hits.begin(), other.hits.end());
    }

    void
    print(std::ostream& out = std::cout) const
    {
      for (const auto& v : hits)
        out << v(0) / units::_cm << " " << v(1) / units::_cm << " "
            << v(2) / units::_cm << std::endl;
    }
  };

  typedef this_result result_type;

  template <typename TrackParameters>
  void
  operator()(const TrackParameters& current,
             const TrackParameters& previous,
             result_type&           result) const
  {
    double r1 = previous.position().perp();
    double z1 = previous.position().z();
    double r2 = current.position().perp();
    double z2 = current.position().z();

    // check for hit in barrel
    for (const auto& t : barrel) {
      double r    = std::get<0>(t);
      double zmin = std::get<1>(t);
      double zmax = std::get<2>(t);
      if ((r1 - r) * (r2 - r) < 0) {
        double s   = (r - r1) / (r2 - r1);
        auto   pos = s * current.position() + (1 - s) * previous.position();
        if ((zmin - pos.z()) * (zmax - pos.z()) < 0) result.hits.push_back(pos);
      }
    }

    // check for hit in endcaps
    for (const auto& t : endcaps) {
      double z    = std::get<0>(t);
      double rmin = std::get<1>(t);
      double rmax = std::get<2>(t);
      if ((z1 - z) * (z2 - z) < 0) {
        double s   = (z - z1) / (z2 - z1);
        auto   pos = s * current.position() + (1 - s) * previous.position();
        if ((rmin - pos.perp()) * (rmax - pos.perp()) < 0)
          result.hits.push_back(pos);
      }
    }
  }

  std::list<std::tuple<float, float, float>> barrel;
  std::list<std::tuple<float, float, float>> endcaps;
};

struct MaxRadius
{
  double maxRadius = 0;

  template <typename TrackParameters>
  bool
  operator()(const TrackParameters& pars, double& stepMax) const
  {
    // adjust maximum step size
    // dR/dS = pT/p -> dS = dR * p / pT = dR * sqrt( 1 + pZ^2 / pT^2)
    const double deltaR       = pars.position().perp() - maxRadius;
    const double pT2          = pars.momentum().perp2();
    const double pZ2          = pars.momentum().z() * pars.momentum().z();
    const double inv_sintheta = sqrt(1 + pZ2 / pT2);

    stepMax = -deltaR * inv_sintheta;

    return (deltaR > 0);
  }
};

int
main()
{
  typedef ConstantBField            BField_type;
  typedef AtlasStepper<BField_type> Stepper_type;
  typedef Propagator<Stepper_type>  Propagator_type;

  BField_type     bField(0, 0, 2 * units::_T);
  Stepper_type    gsl_stepper(std::move(bField));
  Propagator_type propagator(std::move(gsl_stepper));

  // first example: simple propagation of a single track
  {
    typedef ObserverList<DebugObserver> ObsList_type;
    typedef AbortList<MaxRadius>        AbortConditions_type;

    Propagator_type::Options<ObsList_type, AbortConditions_type> options;
    options.max_path_length       = 35 * units::_cm;
    AbortConditions_type& al      = options.stop_conditions;
    al.get<MaxRadius>().maxRadius = 15 * units::_cm;

    std::FILE*    file          = std::fopen("track.txt", "w");
    ObsList_type& ol            = options.observer_list;
    ol.get<DebugObserver>().out = file;
    Vector3D              pos(0, 0, 0);
    Vector3D              mom(1 * units::_GeV, 0, 0);
    CurvilinearParameters pars(nullptr, pos, mom, +1);
    auto                  r = propagator.propagate(pars, options);
    std::cout << "path length = " << r.pathLength / units::_m << std::endl;
    std::fclose(file);
  }

  // second example: use propagator to generate pseudo-hits
  {
    // random numbers for initializing track parameters
    std::random_device                     r;
    std::default_random_engine             e(r());
    std::uniform_real_distribution<double> uniform(0, 1);

    typedef ObserverList<HitSimulator> ObsList_type;

    Propagator_type::Options<ObsList_type> options;
    ObsList_type&                          ol      = options.observer_list;
    HitSimulator&                          hit_sim = ol.get<HitSimulator>();
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
