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

#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Propagator/detail/StandardAbortConditions.hpp"
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
    Vector3D dir = Vector3D(0., 0., 1);
    double   p   = 300. * units::_MeV;

    /// Access method - position
    const Vector3D
    position() const
    {
      return Vector3D(0., 0., 0.);
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
    getField(SteppingState&, const Vector3D&) const
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

  /// @brief mokup of the Options
  struct Options
  {
    bool   loopProtection = true;
    double loopFraction   = 0.5;
  };

  /// @brief mockup of propagtor state
  struct PropagatorState
  {
    /// Contains: stepping state
    SteppingState stepping;
    /// Contains: navigation state
    NavigationState navigation;
    /// Contains: target aborters
    AbortList<PathLimitReached> targetAborters;
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

    PropagatorState pState;
    pState.stepping.dir = Vector3D(cos(phi), sin(phi), 0.);
    pState.stepping.p   = 100. * units::_MeV;

    Stepper pStepper;

    LoopProtection<PathLimitReached> lProtection;
    lProtection(pState, pStepper);

    // pState.targetAborters.get<PathLimitReached>().internalLimit;

    //    // we have a loop logger and a loop aborter
    //    typedef detail::LoopLogger<> LoopLogger;
    //    typedef typename LoopLogger::result_type LoopLoggerResult;
    //    typedef detail::LoopAborter<> LoopAborter;
    //
    //    LoopLogger        loopLogger;
    //    LoopLoggerResult  belowResult;
    //    //LoopLoggerResult  aboveResult;
    //    LoopAborter       loopAborter;
    //
    //    // should not yet be initialized
    //    BOOST_TEST(belowResult.initialized == false);
    //    // first call
    //    loopLogger(below, belowResult);
    //    // should be initialized
    //    BOOST_TEST(belowResult.initialized  == true);
    //    // should not abort yet
    //    loopLogger(below, belowResult);
    //    BOOST_TEST( loopAborter(belowResult, below) == false);
    //    // now move the particle for less than pi
    //    below.stepping.dir = Vector3D(cos(phi+halfPI),sin(phi+halfPI),0.);
    //    // we have not yet met the loop criterium, should not abort yet
    //    loopLogger(below, belowResult);
    //    // test the current delta
    //    BOOST_CHECK_CLOSE(belowResult.currentDelta,halfPI,0.001);
    //    BOOST_TEST( loopAborter(belowResult, below) == false);
    //    // now move the particle for less than pi
    //    below.stepping.dir =
    //    Vector3D(cos(phi+oneAndHalfPI),sin(phi+oneAndHalfPI),0.);
    //    // loop detector should ring a bell
    //    loopLogger(below, belowResult);
    //    // test the current delta
    //    BOOST_CHECK_CLOSE(belowResult.currentDelta,oneAndHalfPI,0.001);
    //    BOOST_TEST( loopAborter(belowResult, below) == true);
    //    // random test
    //    // now move the particle for less than pi
    //    below.stepping.dir = Vector3D(cos(phi+deltaPhi),sin(phi+deltaPhi),0.);
    //    // loop detector should ring a bell
    //    loopLogger(below, belowResult);
    //    // test the current delta
    //    BOOST_CHECK_CLOSE(belowResult.currentDelta,deltaPhi,0.001);
    //    BOOST_TEST( loopAborter(belowResult, below) ==
    //    deltaPhi*deltaPhi>M_PI*M_PI);
  }

}  // namespace Test
}  // namespace Acts