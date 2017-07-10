// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE RungeKuttaEngine Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/EventData/NeutralParameters.hpp"
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/ExtrapolationCell.hpp"
#include "ACTS/Extrapolation/RungeKuttaEngine.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/MagneticField/InterpolatedBFieldMap.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Surfaces/PerigeeSurface.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  BOOST_AUTO_TEST_SUITE(RungeKuttaEngine);

  /// This tests does a propagation through
  /// a magnetic field with 2 Tesla - at nominal incidence
  BOOST_AUTO_TEST_CASE(RungeKuttaEngineTests)
  {

    // use as a test for constant versus interpolated field as well
    // This needs Field.txt from ATLAS
    bool cifTest = false;

    // set up the magnetic field
    // - this is a constant 2 Tesla manetic field
    auto constantField
        = std::make_shared<ConstantBField>(0., 0., 2. * units::_T);

    // RungeKuttaEngine - set up the RungeKuttaEngine - for constant field
    using RungeKuttaEngineCF = Acts::RungeKuttaEngine<ConstantBField>;
    typename RungeKuttaEngineCF::Config rkConfigCF{};
    rkConfigCF.fieldService = constantField;
    auto rkEngineCF         = std::make_shared<RungeKuttaEngineCF>(rkConfigCF);
    rkEngineCF->setLogger(
        getDefaultLogger("RungeKuttaEngineCF", Logging::INFO));

    // RungeKuttaEngine - set up the RungeKuttaEngine - for interpolated field
    using RungeKuttaEngineIF = Acts::RungeKuttaEngine<InterpolatedBFieldMap>;
    std::shared_ptr<RungeKuttaEngineIF> rkEngineIF = nullptr;
    if (cifTest) {
      // - this is a non-constant field
      InterpolatedBFieldMap::Config ncConfig;
      ncConfig.fieldMapFile = "Field.txt";
      auto interpolatedField
          = std::make_shared<InterpolatedBFieldMap>(std::move(ncConfig));
      // set up the engine now
      typename RungeKuttaEngineIF::Config rkConfigIF{};
      rkConfigIF.fieldService = interpolatedField;
      rkEngineIF = std::make_shared<RungeKuttaEngineIF>(rkConfigIF);
      rkEngineIF->setLogger(
          getDefaultLogger("RungeKuttaEngineIF", Logging::INFO));
    }

    // target surface at one meter radius, 10 meters long
    CylinderSurface tSurface(nullptr, 1. * units::_m, 10 * units::_m);

    // the start perigee
    PerigeeSurface pSurface(Vector3D(0., 0., 0.));
    double         d0    = 0. * units::_mm;
    double         z0    = 0. * units::_mm;
    double         phi   = M_PI * 0.25;
    double         theta = M_PI * 0.5;
    double         pT    = 500. * units::_MeV;

    // parameters
    ActsVectorD<5> nparameters;
    nparameters << d0, z0, phi, theta, 1. / pT;

    // some screen output
    std::unique_ptr<Acts::ActsSymMatrixD<5>> cov = nullptr;

    // a straight line propation first
    NeutralBoundParameters nParameters(nullptr, nparameters, pSurface);
    ExtrapolationCell<NeutralParameters> enc(nParameters);

    // go and extrapolate to the
    auto eCode = rkEngineCF->propagate(enc, tSurface);
    /// test for SUCCESS
    int eCodeSuccess = (int)ExtrapolationCode::SuccessDestination;
    BOOST_CHECK_EQUAL(eCode.code, eCodeSuccess);

    // positively charged
    ActsVectorD<5> pcparameters;
    pcparameters << d0, z0, phi, theta, 1. / pT;
    BoundParameters pChargedParameters(nullptr, pcparameters, pSurface);
    ExtrapolationCell<TrackParameters> epcc(pChargedParameters);
    eCode = rkEngineCF->propagate(epcc, tSurface);
    /// test for SUCCESS
    BOOST_CHECK_EQUAL(eCode.code, eCodeSuccess);

    // negatively charged
    ActsVectorD<5> ncparameters;
    ncparameters << d0, z0, phi, theta, -1. / pT;
    BoundParameters nChargedParameters(nullptr, ncparameters, pSurface);
    ExtrapolationCell<TrackParameters> encc(nChargedParameters);
    eCode = rkEngineCF->propagate(encc, tSurface);
    /// test for SUCCESS
    BOOST_CHECK_EQUAL(eCode.code, eCodeSuccess);

    // positively charged - interpolated field
    if (rkEngineIF) {

      ActsVectorD<5> pcparametersIF;
      pcparametersIF << d0, z0, phi, theta, 1. / pT;
      BoundParameters pChargedParametersIF(nullptr, ncparameters, pSurface);
      ExtrapolationCell<TrackParameters> epccIF(pChargedParametersIF);
      eCode = rkEngineIF->propagate(epccIF, tSurface);
      /// test for SUCCESS
      BOOST_CHECK_EQUAL(eCode.code, eCodeSuccess);
    }

    // now check the radii
    if (epcc.endParameters && encc.endParameters) {
      BOOST_TEST_MESSAGE("Testing validity of result.");
      // get the positions
      auto gpc = epcc.endParameters->position();
      auto gpn = encc.endParameters->position();
      // x and y should just be swapped
      BOOST_CHECK_CLOSE(gpc.x(), gpn.y(), 0.001);
      BOOST_CHECK_CLOSE(gpc.y(), gpn.x(), 0.001);
      // and the numbers are
      BOOST_CHECK_CLOSE(gpc.x(), 989.876, 0.001);
      BOOST_CHECK_CLOSE(gpc.y(), 141.935, 0.001);
    } else
      BOOST_TEST_MESSAGE("Propagations failed !");
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test
}  // end of namespace Acts
