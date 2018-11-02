// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE LoopAborter Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using namespace detail;

namespace Test {

  /// @brief mockup of stepping state
  struct SteppingState
  {

    /// Parameters
    Vector3D pos = Vector3D(0., 0., 0.);
    Vector3D dir = Vector3D(0., 0., 1);
    double   p   = 100. * units::_MeV;

    /// Access method - position
    const Vector3D
    position() const
    {
      return pos;
    }

    /// Access method - direction
    const Vector3D
    direction() const
    {
      return dir;
    }

    /// Access method - momentum
    const Vector3D
    momentum() const
    {
      return p * dir;
    }
  };

  /// @brief mockup of stepping state
  struct Stepper
  {

    Vector3D field = Vector3D(0., 0., 2. * units::_T);

    /// Get the field for the stepping, it checks first if the access is still
    /// within the Cell, and updates the cell if necessary.
    ///
    /// @param [in,out] state is the propagation state associated with the track
    ///                 the magnetic field cell is used (and potentially
    ///                 updated)
    /// @param [in] pos is the field position
    Vector3D
    getField(SteppingState& /*unused*/, const Vector3D& /*unused*/) const
    {
      // get the field from the cell
      return field;
    }
  };

  /// @brief mockup of navigation state
  struct NavigationState
  {
    bool navigationBreak = false;
  };

  /// @brief mockup of the Propagator Options
  struct Options
  {

    /// Absolute maximum path length
    double pathLimit      = std::numeric_limits<double>::max();
    bool   loopProtection = true;
    double loopFraction   = 0.5;

    /// Contains: target aborters
    AbortList<PathLimitReached> abortList;
  };

  /// @brief mockup of propagtor state
  struct PropagatorState
  {
    /// Contains: stepping state
    SteppingState stepping;
    /// Contains: navigation state
    NavigationState navigation;
    /// Contains: options
    Options options;
  };

  // This test case checks that no segmentation fault appears
  // - this tests the collection of surfaces
  BOOST_DATA_TEST_CASE(
      loop_aborter_test,
      bdata::random((bdata::seed = 21,
                     bdata::distribution
                     = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::xrange(1),
      phi,
      deltaPhi,
      index)
  {
    (void)index;
    (void)deltaPhi;

    PropagatorState pState;
    pState.stepping.dir = Vector3D(cos(phi), sin(phi), 0.);
    pState.stepping.p   = 100. * units::_MeV;

    Stepper pStepper;

    auto initialLimit
        = pState.options.abortList.get<PathLimitReached>().internalLimit;

    LoopProtection<PathLimitReached> lProtection;
    lProtection(pState, pStepper);

    auto updatedLimit
        = pState.options.abortList.get<PathLimitReached>().internalLimit;
    BOOST_TEST(updatedLimit < initialLimit);
  }

  using BField          = ConstantBField;
  using EigenStepper    = EigenStepper<BField>;
  using EigenPropagator = Propagator<EigenStepper>;

  const int ntests = 100;
  const int skip   = 0;

  // This test case checks that the propagator with loop LoopProtection
  // stops where expected
  BOOST_DATA_TEST_CASE(
      propagator_loop_protection_test,
      bdata::random((bdata::seed = 20,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.5 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 21,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 22,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 23,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    if (index < skip) {
      return;
    }

    double dcharge = -1 + 2 * charge;

    const double Bz = 2. * units::_T;
    BField       bField(0, 0, Bz);

    EigenStepper estepper(bField);

    EigenPropagator epropagator(std::move(estepper));

    // define start parameters
    double                x  = 0;
    double                y  = 0;
    double                z  = 0;
    double                px = pT * cos(phi);
    double                py = pT * sin(phi);
    double                pz = pT / tan(theta);
    double                q  = dcharge;
    Vector3D              pos(x, y, z);
    Vector3D              mom(px, py, pz);
    CurvilinearParameters start(nullptr, pos, mom, q);

    PropagatorOptions<> options;
    const auto&         result = epropagator.propagate(start, options);

    // this test assumes state.options.loopFraction = 0.5
    BOOST_CHECK_CLOSE(px, -result.endParameters->momentum().x(), 1e-2);
    BOOST_CHECK_CLOSE(py, -result.endParameters->momentum().y(), 1e-2);
    BOOST_CHECK_CLOSE(pz, result.endParameters->momentum().z(), 1e-2);
  }

}  // namespace Test
}  // namespace Acts
