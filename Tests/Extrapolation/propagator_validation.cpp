#define BOOST_TEST_MODULE Propagator validation
#include <boost/test/included/unit_test.hpp>

#include <boost/test/data/test_case.hpp>
#include <fstream>
#include <random>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/EigenStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "atlas_propagator_fixture.hpp"
#include "covariance_validation_fixture.hpp"
#include "logfile_erasure_fixture.hpp"

namespace bdata = boost::unit_test::data;
namespace utf   = boost::unit_test;

namespace Acts {

using namespace propagation;

namespace Test {

  BOOST_AUTO_TEST_SUITE(propagator_validation,
                        *utf::fixture<logfile_erasure_fixture>(
                            std::string("gsl_stepper_validation.txt")))

  BOOST_DATA_TEST_CASE_F(
      atlas_propagator_fixture,
      gsl_stepper_validation,
      bdata::random((bdata::seed = 100,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4, 10.)))
          ^ bdata::random((bdata::seed = 100,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., 2 * M_PI)))
          ^ bdata::random((bdata::seed = 100,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::xrange(1000),
      pT,
      phi,
      theta,
      index)
  {
    typedef ConstantBField            BField_type;
    typedef EigenStepper<BField_type> Stepper_type;
    typedef Propagator<Stepper_type>  Propagator_type;

    BField_type     bField(0, 0, 2 * units::_T);
    Stepper_type    gsl_stepper(std::move(bField));
    Propagator_type test_propagator(std::move(gsl_stepper));

    Propagator_type::Options<> options;
    options.max_path_length = 5 * units::_m;

    std::ofstream out("gsl_stepper_validation.txt",
                      std::ios::out | std::ios::app);
    out.precision(6);

    Vector3D pos(0, 0, 0);
    Vector3D mom(pT * cos(phi) * units::_GeV,
                 pT * sin(phi) * units::_GeV,
                 pT * units::_GeV / tan(theta));
    ActsSymMatrixD<5> cov;
    cov << 10, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0.1;

    CurvilinearParameters pars1(
        std::make_unique<ActsSymMatrixD<5>>(cov), pos, mom, +1);
    CurvilinearParameters pars2(
        std::make_unique<ActsSymMatrixD<5>>(cov), pos, mom, +1);

    ExtrapolationCell<TrackParameters> exCell(pars2);
    exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);
    exCell.addConfigurationMode(ExtrapolationMode::CollectJacobians);
    // perform propagation with test propagator
    const auto& r = test_propagator.propagate(pars1, options);
    const auto& p = r.endParameters;
    out << std::fixed << p->position().x() / units::_mm << " "
        << " " << p->position().y() / units::_mm << " "
        << p->position().z() / units::_mm << " "
        << p->momentum().x() / units::_MeV << " "
        << p->momentum().y() / units::_MeV << " "
        << p->momentum().z() / units::_MeV << " " << p->parameters()(0) << " "
        << p->parameters()(1) << " " << p->parameters()(2) << " "
        << p->parameters()(3) << " " << p->parameters()(4) << " "
        << (*p->covariance())(0, 0) << " " << (*p->covariance())(1, 1) << " "
        << (*p->covariance())(2, 2) << " " << (*p->covariance())(3, 3) << " "
        << (*p->covariance())(4, 4) << " ";

    // perform propagation with ATLAS Runge-Kutta propagator
    exCell.leadParameters     = &pars2;
    exCell.lastLeadParameters = 0;
    propagator->propagate(exCell, *surface);
    const auto* val = exCell.leadParameters;
    out << std::fixed << val->position().x() / units::_mm << " "
        << " " << val->position().y() / units::_mm << " "
        << val->position().z() / units::_mm << " "
        << val->momentum().x() / units::_MeV << " "
        << val->momentum().y() / units::_MeV << " "
        << val->momentum().z() / units::_MeV << " " << val->parameters()(0)
        << " " << val->parameters()(1) << " " << val->parameters()(2) << " "
        << val->parameters()(3) << " " << val->parameters()(4) << " "
        << (*val->covariance())(0, 0) << " " << (*val->covariance())(1, 1)
        << " " << (*val->covariance())(2, 2) << " "
        << (*val->covariance())(3, 3) << " " << (*val->covariance())(4, 4);

    out << " " << r.pathLength << " " << exCell.pathLength << std::endl;
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // namespace Test
}  // namespace Acts
