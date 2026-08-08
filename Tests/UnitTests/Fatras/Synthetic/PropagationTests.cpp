// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/detail/Helix.hpp"
#include "ActsFatras/Synthetic/detail/Propagation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <vector>

using namespace ActsFatras::Synthetic;
using namespace ActsFatras::Synthetic::detail;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

constexpr float kBFieldZ = 2_T;

/// Enough for a track to cross both kinds of surface.
DetectorLayout makeTestLayout() {
  DetectorLayoutBuilder builder;
  builder.addPassiveCylinder(25.f)
      .addCylinder(34.f, 250.f, 1)
      .addCylinder(99.f, 250.f, 1)
      .addCylinder(160.f, 250.f, 1);
  for (const SurfaceSide side :
       {SurfaceSide::Positive, SurfaceSide::Negative}) {
    for (const float absZ : {400.f, 600.f, 800.f}) {
      builder.addDisc(side, absZ, 30.f, 200.f, 1);
    }
  }
  return builder.build();
}

/// Launch a secondary from a crossing of its parent, the way the generator
/// does.
Helix makeSecondary(const Helix& parent, const SurfaceCrossing& crossing) {
  const float direction = parent.phi0 - parent.charge * crossing.gamma;
  return makeHelixFromPoint(crossing.r, crossing.phi, crossing.z, direction,
                            1.2f * parent.cotTheta, 0.4_GeV, -parent.charge,
                            kBFieldZ);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticPropagationSuite)

/// Every crossing is a point on the helix at the angle it reports, inside the
/// layer it resolved to.
BOOST_AUTO_TEST_CASE(CrossingsLieOnTheHelix) {
  constexpr float pi = std::numbers::pi_v<float>;
  const DetectorLayout layout = makeTestLayout();

  std::vector<SurfaceCrossing> crossings;
  std::size_t inbound = 0;
  for (const float cotTheta : {0.2f, 1.5f, 3.f}) {
    for (const float charge : {-1.f, +1.f}) {
      const Helix helix =
          makeHelix(0.1f, 0.3f, 5.f, cotTheta, 2_GeV, charge, kBFieldZ);
      crossings.clear();
      // more than a half turn, so the inbound legs are exercised too
      propagateHelix(layout, helix, crossings, kNoSurface, 2.f * pi);
      BOOST_REQUIRE(!crossings.empty());

      std::uint32_t previous = 0;
      for (const SurfaceCrossing& crossing : crossings) {
        BOOST_REQUIRE_LT(crossing.surface, layout.surfaces.size());

        // the position is what the helix has at that turning angle
        BOOST_CHECK_CLOSE(crossing.r, helixRadiusAtGamma(helix, crossing.gamma),
                          1e-2);
        const float z =
            helix.z0 + helix.cotTheta * helix.radius * crossing.gamma;
        BOOST_CHECK_SMALL(crossing.z - z, 1e-2f);
        const float phi = helixPhiAtRadius(helix, crossing.r);
        const float expected =
            helixIsInbound(crossing.gamma) ? helixPhiMirrored(helix, phi) : phi;
        BOOST_CHECK_SMALL(std::remainder(crossing.phi - expected, 2.f * pi),
                          1e-3f);
        inbound += helixIsInbound(crossing.gamma) ? 1 : 0;

        // a sensitive crossing resolves to a layer that contains it
        const DetectorSurface& surface = layout.surfaces[crossing.surface];
        if (crossing.sensitive()) {
          BOOST_REQUIRE_LT(crossing.layer, layout.layers.size());
          const DetectorLayer& layer = layout.layers[crossing.layer];
          const float along =
              layer.shape == SurfaceShape::Cylinder ? crossing.z : crossing.r;
          BOOST_CHECK_GE(along, layer.minBound);
          BOOST_CHECK_LE(along, layer.maxBound);
        } else {
          BOOST_CHECK(surface.layers.empty());
        }

        // and the crossings come in surface order
        BOOST_CHECK_GE(crossing.surface, previous);
        previous = crossing.surface;
      }
    }
  }
  // the mirrored azimuth above was actually reached
  BOOST_CHECK_GT(inbound, 0u);
}

/// A track produced on a surface does not cross it: its turning angle there
/// *is* `Helix::minGamma`, so only naming the surface can exclude it. What goes
/// is the departure alone -- a secondary soft enough to curl back through the
/// surface that made it still measures it on the way in.
BOOST_AUTO_TEST_CASE(ProductionSurfaceSkipsOnlyTheDeparture) {
  constexpr float pi = std::numbers::pi_v<float>;
  const DetectorLayout layout = makeTestLayout();

  std::vector<SurfaceCrossing> crossings;
  std::vector<SurfaceCrossing> all;
  std::vector<SurfaceCrossing> rest;
  std::size_t discs = 0;
  std::size_t cylinders = 0;
  std::size_t returned = 0;

  // A scan rather than a handful of tracks: whether the production disc came
  // back out was decided by rounding, so only some configurations showed it.
  for (std::size_t i = 0; i < 60; ++i) {
    const float cotTheta = 0.2f + 0.1f * static_cast<float>(i);
    for (const float charge : {-1.f, +1.f}) {
      const Helix parent =
          makeHelix(0.2f, -0.6f, 12.f, cotTheta, 1.5_GeV, charge, kBFieldZ);
      crossings.clear();
      propagateHelix(layout, parent, crossings);

      for (const SurfaceCrossing& crossing : crossings) {
        // soft, and followed for turns, so that it comes back through the
        // surface it was made on
        const Helix secondary = makeSecondary(parent, crossing);
        const float maxGamma = secondary.minGamma + 4.f * pi;
        all.clear();
        propagateHelix(layout, secondary, all, kNoSurface, maxGamma);
        rest.clear();
        propagateHelix(layout, secondary, rest, crossing.surface, maxGamma);

        // What goes is the crossing at the production angle, and nothing else.
        // Not an exact comparison against `minGamma`: on a disc the departure
        // is solved for rather than carried over, and rounds either side of it.
        const auto departure = [&](const SurfaceCrossing& other) {
          return other.surface == crossing.surface &&
                 std::abs(other.gamma - secondary.minGamma) < 1e-4f;
        };
        const auto departures =
            static_cast<std::size_t>(std::ranges::count_if(all, departure));
        BOOST_REQUIRE_EQUAL(rest.size(), all.size() - departures);
        BOOST_CHECK(std::ranges::none_of(rest, departure));
        returned += std::ranges::count_if(rest, [&](const auto& other) {
          return other.surface == crossing.surface;
        });

        if (layout.surfaces[crossing.surface].shape == SurfaceShape::Disc) {
          ++discs;
        } else {
          ++cylinders;
        }
      }
    }
  }

  // both kinds of production surface were exercised, and the re-crossing the
  // check above has to keep really does happen
  BOOST_CHECK_GT(discs, 0u);
  BOOST_CHECK_GT(cylinders, 0u);
  BOOST_CHECK_GT(returned, 0u);
}

/// A crossing that lands off every ring of a disc still meets the support.
BOOST_AUTO_TEST_CASE(MaterialIsMetOffTheRings) {
  const Acts::MaterialSlab slab =
      materialSlab(95.7f, 465.2f, 28.03f, 14.f, 8.28e-8f, 1.4f);
  DetectorLayoutBuilder builder;
  const std::array<RingBounds, 2> rings{RingBounds{40.f, 60.f},
                                        RingBounds{100.f, 140.f}};
  builder.addDisc(SurfaceSide::Positive, 600.f, rings);
  DetectorLayout layout = builder.build();
  // material over the whole disc, rings and the gap between them alike
  layout.surfaces.front().material =
      SurfaceMaterial{BandComposition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f},
                      {30.f, 150.f},
                      {0.f, 95.7f},
                      {0.f, 465.2f}};

  // straight enough to cross the disc between the two rings
  const Helix helix =
      makeHelix(0.f, 0.f, 0.f, 600.f / 80.f, 100_GeV, -1.f, kBFieldZ);
  std::vector<SurfaceCrossing> crossings;
  propagateHelix(layout, helix, crossings);

  BOOST_REQUIRE_EQUAL(crossings.size(), 1u);
  BOOST_CHECK(!crossings.front().sensitive());
  BOOST_CHECK_CLOSE(crossings.front().r, 80.f, 1e-2);
  BOOST_CHECK_CLOSE(crossings.front().material->thicknessInX0(),
                    slab.thicknessInX0(), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
