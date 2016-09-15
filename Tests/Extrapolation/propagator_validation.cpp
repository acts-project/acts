#define BOOST_TEST_MODULE Propagator validation
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/AbortConditions.hpp"
#include "ACTS/Extrapolation/AbortList.hpp"
#include "ACTS/Extrapolation/GSLStepper.hpp"
#include "ACTS/Extrapolation/ObserverList.hpp"
#include "ACTS/Extrapolation/Observers.hpp"
#include "ACTS/Extrapolation/Propagator.hpp"
#include "ACTS/MagneticField/ConstantFieldSvc.hpp"
#include "ACTS/Utilities/Units.hpp"
#include "atlas_propagator_fixture.hpp"

namespace Acts {

namespace Test {

  BOOST_FIXTURE_TEST_CASE(gsl_stepper_validation, atlas_propagator_fixture)
  {
    typedef ConstantFieldSvc         BField_type;
    typedef GSLStepper<BField_type>  Stepper_type;
    typedef Propagator<Stepper_type> Propagator_type;

    BField_type::Config c;
    c.field = {0, 0, 2 * units::_T};
    BField_type     bField(std::move(c));
    Stepper_type    gsl_stepper(std::move(bField));
    Propagator_type test_propagator(std::move(gsl_stepper));

    typedef ObserverList<PathLengthObserver> ObsList_type;
    typedef AbortList<MaxPathLength>         AbortConditions_type;

    Propagator_type::Options<ObsList_type, AbortConditions_type> options;
    AbortConditions_type& al              = options.stop_conditions;
    al.get<MaxPathLength>().maxPathLength = 35 * units::_cm;

    std::ofstream out("gsl_stepper_validation.txt");
    out.precision(4);

    Vector3D              pos(0, 0, 0);
    Vector3D              mom(1 * units::_GeV, 0, 0);
    CurvilinearParameters pars1(nullptr, pos, mom, +1);
    CurvilinearParameters pars2(nullptr, pos, mom / units::_MeV, +1);

    ExtrapolationCell<TrackParameters> exCell(pars2);
    exCell.addConfigurationMode(ExtrapolationMode::StopWithPathLimit);
    for (unsigned int i = 0; i < 100; ++i) {
      // perform propagation with test propagator
      const auto& p = test_propagator.propagate(pars1, options).endParameters;
      out << std::fixed << p.position().x() / units::_mm << " "
          << " " << p.position().y() / units::_mm << " "
          << p.position().z() / units::_mm << " "
          << p.momentum().x() / units::_GeV << " "
          << p.momentum().y() / units::_GeV << " "
          << p.momentum().z() / units::_GeV << " ";

      // perform propagation with ATLAS Runge-Kutta propagator
      exCell.leadParameters     = &pars2;
      exCell.lastLeadParameters = 0;
      propagator->propagate(exCell, *surface);
      const auto* val = exCell.leadParameters;
      out << std::fixed << val->position().x() / units::_mm << " "
          << " " << val->position().y() / units::_mm << " "
          << val->position().z() / units::_mm << " "
          << val->momentum().x() * units::_MeV << " "
          << val->momentum().y() * units::_MeV << " "
          << val->momentum().z() * units::_MeV << std::endl;
    }
  }
}  // namespace Test
}  // namespace Acts
