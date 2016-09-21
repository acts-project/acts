#define BOOST_TEST_MODULE Propagator validation
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <fstream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/AtlasStepper.hpp"
#include "ACTS/Extrapolation/GSLStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Observers.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "atlas_propagator_fixture.hpp"
#include "logfile_erasure_fixture.hpp"

namespace bdata = boost::unit_test::data;
namespace utf   = boost::unit_test;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_SUITE(propagator_validation,
                        *utf::fixture<logfile_erasure_fixture>(
                            std::string("gsl_stepper_validation.txt")))

  BOOST_DATA_TEST_CASE_F(atlas_propagator_fixture,
                         gsl_stepper_validation,
                         bdata::random(0.4, 10.) ^ bdata::random(0., 2 * M_PI)
                             ^ bdata::random(0., M_PI)
                             ^ bdata::xrange(10000),
                         pT,
                         phi,
                         theta,
                         index)
  {
    typedef ConstantFieldSvc          BField_type;
    typedef AtlasStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    BField_type::Config c;
    c.field = {0, 0, 2 * units::_T};
    BField_type     bField(std::move(c));
    Stepper_type    gsl_stepper(std::move(bField));
    Propagator_type test_propagator(std::move(gsl_stepper));

    typedef ObserverList<PathLengthObserver> ObsList_type;
    typedef AbortList<MaxPathLength>         AbortConditions_type;

    Propagator_type::Options<ObsList_type, AbortConditions_type> options;
    AbortConditions_type& al              = options.stop_conditions;
    al.get<MaxPathLength>().maxPathLength = 5 * units::_m;

    std::ofstream out("gsl_stepper_validation.txt",
                      std::ios::out | std::ios::app);
    out.precision(6);

    Vector3D pos(0, 0, 0);
    Vector3D mom(pT * cos(phi) * units::_GeV,
                 pT * sin(phi) * units::_GeV,
                 pT * units::_GeV / tan(theta));
    CurvilinearParameters pars1(nullptr, pos, mom, +1);
    CurvilinearParameters pars2(nullptr, pos, mom / units::_MeV, +1);

    ExtrapolationCell<TrackParameters> exCell(pars2);
    exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);
    // perform propagation with test propagator
    const auto& r = test_propagator.propagate(pars1, options);
    const auto& p = r.endParameters;
    out << std::fixed << p.position().x() / units::_mm << " "
        << " " << p.position().y() / units::_mm << " "
        << p.position().z() / units::_mm << " "
        << p.momentum().x() / units::_MeV << " "
        << p.momentum().y() / units::_MeV << " "
        << p.momentum().z() / units::_MeV << " " << p.parameters()(0) << " "
        << p.parameters()(1) << " " << p.parameters()(2) << " "
        << p.parameters()(3) << " " << p.parameters()(4) << " ";

    // perform propagation with ATLAS Runge-Kutta propagator
    exCell.leadParameters     = &pars2;
    exCell.lastLeadParameters = 0;
    propagator->propagate(exCell, *surface);
    const auto* val = exCell.leadParameters;
    out << std::fixed << val->position().x() / units::_mm << " "
        << " " << val->position().y() / units::_mm << " "
        << val->position().z() / units::_mm << " " << val->momentum().x() << " "
        << val->momentum().y() << " " << val->momentum().z() << " "
        << val->parameters()(0) << " " << val->parameters()(1) << " "
        << val->parameters()(2) << " " << val->parameters()(3) << " "
        << val->parameters()(4) * 1000;  //<< std::endl;

    out << " " << r.get<PathLengthObserver::result_type>().pathLength << " "
        << exCell.pathLength << std::endl;
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace Test
}  // namespace Acts
