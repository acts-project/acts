// This file is part of the ACTS project.
//
// Copyright (C) 2016-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE Propagator Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/MagneticField/ConstantBField.hpp"
#include "ACTS/Propagator/EigenStepper.hpp"
#include "ACTS/Propagator/ObserverList.hpp"
#include "ACTS/Propagator/Propagator.hpp"
#include "ACTS/Propagator/detail/standard_abort_conditions.hpp"
#include "ACTS/Surfaces/CylinderSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

using namespace propagation;

namespace Test {

  /// An observer that measures the perpendicular distance
  struct PerpendicularMeasure
  {

    /// Simple result struct to be returned
    struct this_result
    {
      double distance = std::numeric_limits<double>::max();
    };

    typedef this_result result_type;

    PerpendicularMeasure() {}

    template <typename input_t>
    void
    operator()(const input_t& cache, result_type& result) const
    {
      result.distance = cache.pos.perp();
    }

    template <typename input_t>
    void
    operator()(const input_t& cache) const
    {
      (void)cache;
    }
  };

  // Global definitions
  // The path limit abort
  typedef detail::path_limit_reached path_limit;

  typedef ConstantBField                BField_type;
  typedef EigenStepper<BField_type>     EigenStepper_type;
  typedef Propagator<EigenStepper_type> EigenPropagator_type;

  const double         Bz = 2. * units::_T;
  BField_type          bField(0, 0, Bz);
  EigenStepper_type    estepper(bField);
  EigenPropagator_type epropagator(std::move(estepper));

  CylinderSurface cSurface(nullptr, 150., 1000. * units::_mm);

  // This tests the Options
  BOOST_AUTO_TEST_CASE(PropgatorOptions_)
  {

    typedef typename Propagator<EigenStepper_type>::template Options<>
                      null_options_type;
    null_options_type null_options;
    // todo write null options test

    typedef ObserverList<PerpendicularMeasure> ObsList_type;
    typedef AbortList<>                        AbortConditions_type;

    typedef typename Propagator<EigenStepper_type>::
        template Options<ObsList_type, AbortConditions_type>
            options_type;

    options_type options;

  }

  BOOST_DATA_TEST_CASE(
      constant_bfield_cylinder_intersection_observer_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(100),
      pT,
      phi,
      theta,
      charge,
      index)
  {
    double dcharge = -1 + 2 * charge;
    // constant field propagation atlas stepper
    (void)index;

    typedef ObserverList<PerpendicularMeasure> ObsList_type;
    typedef AbortList<>                        AbortConditions_type;

    //   // setup propagation options
    typename EigenPropagator_type::template Options<ObsList_type,
                                                    AbortConditions_type>
        options;

    options.max_path_length = 20 * units::_m;
    options.max_step_size   = 1 * units::_cm;

    typedef typename PerpendicularMeasure::result_type pm_result;

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
    // propagate to the cylinder surface
    const auto& result = epropagator.propagate(start, cSurface, options);
    auto&       pmr    = result.get<pm_result>();

    // Test the end position
    BOOST_TEST(std::abs(result.endParameters->position().perp() - 150.)
               < 10e-4);
    BOOST_TEST(std::abs(pmr.distance - 150.) < 10e-4);
  }

}  // namespace Test
}  // namespace Acts
