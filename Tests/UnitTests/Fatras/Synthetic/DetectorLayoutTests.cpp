// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/detail/Propagation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace ActsFatras::Synthetic;
using ActsFatras::Synthetic::detail::findLayer;

namespace ActsTests {

namespace {

/// The invariants every layout has to satisfy: the layers of a surface increase
/// along its free coordinate without overlapping, agree with it on shape and
/// reference coordinate, and belong to no other surface.
void checkLayoutInvariants(const DetectorLayout& layout) {
  std::set<std::uint32_t> seen;
  for (const DetectorSurface& surface : layout.surfaces) {
    // whatever the surface holds is inside its extent, or propagation would
    // reject a crossing that lands on it
    BOOST_CHECK_LE(surface.minBound, surface.maxBound);
    if (surface.layers.empty()) {
      // passive, nothing to tile; it was declared with its extent
      continue;
    }
    for (std::size_t k = 0; k < surface.material.bounds.size(); ++k) {
      if (surface.material.bands[k].isVacuum()) {
        continue;
      }
      BOOST_CHECK_LE(surface.minBound,
                     k == 0 ? 0.f : surface.material.bounds[k - 1]);
      BOOST_CHECK_GE(surface.maxBound, surface.material.bounds[k]);
    }
    for (const std::uint32_t index : surface.layers) {
      const DetectorLayer& layer = layout.layers[index];
      // a cylinder is bounded in |z|, and its layers are given in signed z
      const float lo =
          surface.shape == SurfaceShape::Cylinder
              ? std::min(std::abs(layer.minBound), std::abs(layer.maxBound))
              : layer.minBound;
      const float hi =
          surface.shape == SurfaceShape::Cylinder
              ? std::max(std::abs(layer.minBound), std::abs(layer.maxBound))
              : layer.maxBound;
      BOOST_CHECK_LE(surface.minBound, lo);
      BOOST_CHECK_GE(surface.maxBound, hi);
    }

    float previousMax = std::numeric_limits<float>::quiet_NaN();
    std::uint32_t expectedModule = 0;
    for (const std::uint32_t index : surface.layers) {
      BOOST_REQUIRE_LT(index, layout.layers.size());
      // no layer is claimed twice
      BOOST_CHECK(seen.insert(index).second);

      const DetectorLayer& layer = layout.layers[index];
      BOOST_CHECK(layer.shape == surface.shape);
      BOOST_CHECK_EQUAL(layer.refCoord, surface.refCoord);
      BOOST_CHECK_LT(layer.minBound, layer.maxBound);
      BOOST_CHECK_EQUAL(layer.moduleIndex, expectedModule);
      ++expectedModule;

      if (!std::isnan(previousMax)) {
        if (surface.shape == SurfaceShape::Cylinder) {
          // the previous module ends exactly where this one starts
          BOOST_CHECK_EQUAL(layer.minBound, previousMax);
        } else {
          BOOST_CHECK_GE(layer.minBound, previousMax);
        }
      }
      previousMax = layer.maxBound;
    }

    // a cylinder is symmetric in z, a disc sits on one side
    if (surface.shape == SurfaceShape::Cylinder) {
      BOOST_CHECK(layout.layers[surface.layers.front()].side ==
                  SurfaceSide::Barrel);
    } else {
      const SurfaceSide side = layout.layers[surface.layers.front()].side;
      BOOST_CHECK(side != SurfaceSide::Barrel);
      BOOST_CHECK_EQUAL(side == SurfaceSide::Positive, surface.refCoord > 0.f);
    }
  }

  // every layer belongs to some surface
  BOOST_CHECK_EQUAL(seen.size(), layout.layers.size());
}

}  // namespace

namespace {
/// How many surfaces a description's services build: a disc is mirrored onto
/// both sides, a cylinder is not.
std::size_t passiveSurfaceCount(const BarrelEndcapDescription& description) {
  std::size_t count = 0;
  for (const PassiveSurfaceDescription& surface : description.passiveSurfaces) {
    count += surface.shape == SurfaceShape::Cylinder ? 1u : 2u;
  }
  return count;
}
}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticDetectorLayoutSuite)

/// The modules of a surface tile it exactly, and are numbered in call order:
/// cylinders across the whole layout, discs within one endcap.
BOOST_AUTO_TEST_CASE(BuilderTilesAndNumbers) {
  DetectorLayoutBuilder builder;
  builder.addPassiveCylinder(25.f)
      .addCylinder(34.f, 250.f, 4)
      .addCylinder(99.f, 250.f, 1)
      .addDisc(SurfaceSide::Positive, 400.f, 30.f, 350.f, 3)
      .addDisc(SurfaceSide::Positive, 600.f, 30.f, 350.f, 1)
      .addDisc(SurfaceSide::Negative, 400.f, 30.f, 350.f, 3);
  const DetectorLayout layout = builder.build();

  BOOST_CHECK_EQUAL(layout.surfaces.size(), 6u);
  BOOST_CHECK_EQUAL(layout.layers.size(), 4u + 1u + 3u + 1u + 3u);
  checkLayoutInvariants(layout);

  // the split cylinder covers exactly its half-length
  const DetectorSurface& split = layout.surfaces[1];
  BOOST_CHECK_EQUAL(layout.layers[split.layers.front()].minBound, -250.f);
  BOOST_CHECK_EQUAL(layout.layers[split.layers.back()].maxBound, 250.f);
  // and the disc its radial range
  const DetectorSurface& disc = layout.surfaces[3];
  BOOST_CHECK_EQUAL(layout.layers[disc.layers.front()].minBound, 30.f);
  BOOST_CHECK_EQUAL(layout.layers[disc.layers.back()].maxBound, 350.f);

  // cylinders are numbered in call order
  BOOST_CHECK_EQUAL(layout.layers[split.layers.front()].layer, 0u);
  BOOST_CHECK_EQUAL(layout.layers[layout.surfaces[2].layers.front()].layer, 1u);
  // discs are numbered per side, so the first disc of each side is index zero
  BOOST_CHECK_EQUAL(layout.layers[disc.layers.front()].layer, 0u);
  BOOST_CHECK(layout.layers[disc.layers.front()].side == SurfaceSide::Positive);
  BOOST_CHECK_EQUAL(layout.layers[layout.surfaces[4].layers.front()].layer, 1u);
  const DetectorLayer& negative = layout.layers[layout.surfaces[5].layers[0]];
  BOOST_CHECK_EQUAL(negative.layer, 0u);
  BOOST_CHECK(negative.side == SurfaceSide::Negative);
}

BOOST_AUTO_TEST_CASE(BuilderResets) {
  DetectorLayoutBuilder builder;
  builder.addCylinder(34.f, 100.f, 1);
  const DetectorLayout first = builder.build();
  BOOST_CHECK_EQUAL(first.surfaces.size(), 1u);

  // build() hands the layout over and leaves the builder empty, including its
  // layer counters
  builder.addCylinder(70.f, 100.f, 1);
  const DetectorLayout second = builder.build();
  BOOST_CHECK_EQUAL(second.surfaces.size(), 1u);
  BOOST_CHECK_EQUAL(second.layers.size(), 1u);
  BOOST_CHECK_EQUAL(second.layers[0].layer, 0u);
}

BOOST_AUTO_TEST_CASE(BuilderRejectsBadInput) {
  DetectorLayoutBuilder builder;
  BOOST_CHECK_THROW(builder.addPassiveCylinder(0.f), std::invalid_argument);
  // a disc is on an endcap, and the enum keeps everything but the barrel and a
  // value cast in from outside its range out of the argument
  BOOST_CHECK_THROW(builder.addPassiveDisc(SurfaceSide::Barrel, 400.f),
                    std::invalid_argument);
  BOOST_CHECK_THROW(builder.addPassiveDisc(static_cast<SurfaceSide>(+2), 400.f),
                    std::invalid_argument);
  BOOST_CHECK_THROW(builder.addCylinder(34.f, 250.f, 0), std::invalid_argument);
  BOOST_CHECK_THROW(builder.addCylinder(34.f, -1.f, 1), std::invalid_argument);
  BOOST_CHECK_THROW(
      builder.addDisc(SurfaceSide::Positive, 400.f, 350.f, 30.f, 1),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      builder.addDisc(SurfaceSide::Positive, 400.f, 30.f, 350.f, 0),
      std::invalid_argument);
  BOOST_CHECK_THROW(
      builder.addDisc(SurfaceSide::Positive, -400.f, 30.f, 350.f, 1),
      std::invalid_argument);
  // explicit rings: none at all, running backwards, or overlapping
  const std::array backwards{RingBounds{50.f, 30.f}};
  const std::array overlapping{RingBounds{30.f, 60.f}, RingBounds{50.f, 80.f}};
  BOOST_CHECK_THROW(builder.addDisc(SurfaceSide::Positive, 400.f, {}),
                    std::invalid_argument);
  BOOST_CHECK_THROW(builder.addDisc(SurfaceSide::Positive, 400.f, backwards),
                    std::invalid_argument);
  BOOST_CHECK_THROW(builder.addDisc(SurfaceSide::Positive, 400.f, overlapping),
                    std::invalid_argument);
}

/// A crossing that lands in the radial gap between two rings finds no layer.
BOOST_AUTO_TEST_CASE(BuilderDiscWithGaps) {
  DetectorLayoutBuilder builder;
  const std::array rings{RingBounds{154.2f, 194.9f}, RingBounds{214.3f, 255.f},
                         RingBounds{274.3f, 315.1f}};
  builder.addDisc(SurfaceSide::Positive, 1145.5f, rings);
  const DetectorLayout layout = builder.build();

  BOOST_CHECK_EQUAL(layout.surfaces.size(), 1u);
  BOOST_CHECK_EQUAL(layout.layers.size(), 3u);
  checkLayoutInvariants(layout);

  const DetectorSurface& disc = layout.surfaces.front();
  // on a ring
  const std::optional<std::uint32_t> innermost =
      findLayer(layout, disc, 180.f, 1145.5f);
  BOOST_REQUIRE(innermost.has_value());
  BOOST_CHECK_EQUAL(layout.layers[*innermost].moduleIndex, 0u);
  const std::optional<std::uint32_t> outermost =
      findLayer(layout, disc, 300.f, 1145.5f);
  BOOST_REQUIRE(outermost.has_value());
  BOOST_CHECK_EQUAL(layout.layers[*outermost].moduleIndex, 2u);
  // inside the beam hole, in the two gaps, and beyond the outermost ring
  BOOST_CHECK(!findLayer(layout, disc, 100.f, 1145.5f).has_value());
  BOOST_CHECK(!findLayer(layout, disc, 205.f, 1145.5f).has_value());
  BOOST_CHECK(!findLayer(layout, disc, 265.f, 1145.5f).has_value());
  BOOST_CHECK(!findLayer(layout, disc, 400.f, 1145.5f).has_value());
}

/// Every surface takes the material and the overlap its description names. The
/// composition is part of it: the nuclear length is what the hadronic half of
/// the yield counts in, so a carbon support and a silicon sensor cannot share
/// one number.
BOOST_AUTO_TEST_CASE(DescriptionReachesTheSurfaces) {
  const Acts::MaterialSlab silicon =
      materialSlab(95.7f, 465.2f, 28.03f, 14.f, 8.28e-8f, 1.4f);
  const Acts::MaterialSlab carbon =
      materialSlab(284.7f, 572.0f, 12.011f, 6.f, 1.25e-4f, 2.f);
  const Acts::MaterialSlab beryllium =
      materialSlab(352.8f, 421.0f, 9.012f, 4.f, 2.05e-4f, 0.8f);

  BarrelEndcapDescription description;
  description.beamPipeRadius = 25.f;
  description.beamPipeMaterial = SurfaceMaterial{beryllium};
  description.barrelRadii = {34.f, 99.f};
  description.barrelHalfLengthsZ = {244.f, 244.f};
  description.barrelMaterials = {SurfaceMaterial{silicon},
                                 SurfaceMaterial{carbon}};
  description.barrelOverlapProbabilities = {0.02f, 0.15f};
  description.barrelOverlapOffset = 7.9f;
  description.discs = {{600.f, {{40.f, 190.f}}, SurfaceMaterial{carbon}}};
  description.discOverlapProbability = 0.15f;
  description.discOverlapOffset = 5.f;
  description.passiveSurfaces = {
      {SurfaceShape::Cylinder, 60.f, 0.f, 0.f, SurfaceMaterial{carbon}}};

  const DetectorLayout layout = makeLayout(description);
  for (const DetectorSurface& surface : layout.surfaces) {
    const Acts::MaterialSlab& want = surface.refCoord == 25.f   ? beryllium
                                     : surface.refCoord == 34.f ? silicon
                                                                : carbon;
    BOOST_CHECK_CLOSE(surface.material.average().thicknessInX0(),
                      want.thicknessInX0(), 1e-3);
    BOOST_CHECK_CLOSE(surface.material.average().thicknessInL0(),
                      want.thicknessInL0(), 1e-3);

    if (surface.shape == SurfaceShape::Disc) {
      BOOST_CHECK_CLOSE(surface.overlapProbability, 0.15f, 1e-3);
      BOOST_CHECK_CLOSE(surface.overlapOffset, 5.f, 1e-3);
    } else if (surface.layers.empty()) {
      // a beam pipe and a service measure nothing, so they overlap nothing
      BOOST_CHECK_EQUAL(surface.overlapProbability, 0.f);
    } else {
      BOOST_CHECK_CLOSE(surface.overlapProbability,
                        surface.refCoord == 34.f ? 0.02f : 0.15f, 1e-3);
      BOOST_CHECK_CLOSE(surface.overlapOffset, 7.9f, 1e-3);
    }
  }
  // and the two lengths really do disagree, or the check above is vacuous:
  // carbon is worth twice the nuclear length per radiation length silicon is
  BOOST_CHECK_GT(carbon.thicknessInL0() / carbon.thicknessInX0(),
                 1.5f * silicon.thicknessInL0() / silicon.thicknessInX0());
}

/// A banded surface answers with the band the position falls in and with
/// nothing past the last, which is what tells a ring from the gap beside it.
BOOST_AUTO_TEST_CASE(BandedMaterialVariesAlongTheSurface) {
  const Acts::MaterialSlab silicon =
      materialSlab(95.7f, 465.2f, 28.03f, 14.f, 8.28e-8f, 1.4f);
  const SurfaceMaterial uniform{silicon};
  BOOST_CHECK_CLOSE(uniform.at(0.f).thicknessInX0(), silicon.thicknessInX0(),
                    1e-3);
  BOOST_CHECK_CLOSE(uniform.at(1e6f).thicknessInX0(), silicon.thicknessInX0(),
                    1e-3);

  // an endcap disc: a beam hole, two rings of silicon and the carbon support
  // between them, holding a tenth of the material at half the ratio
  const SurfaceMaterial banded{
      BandComposition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f},
      {30.f, 50.f, 80.f, 120.f},
      {0.f, 95.7f, 957.f, 95.7f},
      {0.f, 465.2f, 1914.f, 465.2f}};
  BOOST_CHECK(banded.at(10.f).isVacuum());
  BOOST_CHECK_CLOSE(banded.at(40.f).thicknessInX0(), silicon.thicknessInX0(),
                    1e-3);
  BOOST_CHECK_CLOSE(banded.at(60.f).thicknessInX0(),
                    0.1f * silicon.thicknessInX0(), 1e-3);
  // and it really is made of something else, which is the whole point
  BOOST_CHECK_CLOSE(
      banded.at(60.f).material().L0() / banded.at(60.f).material().X0(), 2.,
      1e-3);
  BOOST_CHECK_CLOSE(
      banded.at(40.f).material().L0() / banded.at(40.f).material().X0(),
      465.2 / 95.7, 1e-3);
  // every band is quoted at the one thickness; the material carries the rest
  BOOST_CHECK_CLOSE(banded.at(60.f).thickness(), banded.at(40.f).thickness(),
                    1e-3);
  // past the outermost band the surface has ended
  BOOST_CHECK(banded.at(130.f).isVacuum());

  // the two lengths per band are a table encoding; giving the bands outright is
  // the same surface, and then the average is taken over them
  const SurfaceMaterial fromSlabs{
      {30.f, 50.f, 80.f, 120.f},
      {Acts::MaterialSlab::Vacuum(silicon.thickness()), silicon,
       materialSlab(957.f, 1914.f, 28.03f, 14.f, 8.28e-9f, 1.4f), silicon}};
  for (const float along : {10.f, 40.f, 60.f, 100.f, 130.f}) {
    BOOST_CHECK_CLOSE(fromSlabs.at(along).thicknessInX0(),
                      banded.at(along).thicknessInX0(), 1e-3);
    BOOST_CHECK_CLOSE(fromSlabs.at(along).thicknessInL0(),
                      banded.at(along).thicknessInL0(), 1e-3);
  }
  BOOST_CHECK_CLOSE(fromSlabs.average().thicknessInX0(),
                    0.25f * (2.1f * silicon.thicknessInX0()), 1e-3);

  // and vacuum stays vacuum however it is asked
  const SurfaceMaterial none;
  BOOST_CHECK(none.at(0.f).isVacuum());
  BOOST_CHECK(none.at(100.f).isVacuum());
}

/// A surface knows where it holds anything, which is what lets propagation
/// answer a crossing of its plane without looking at its layers or its bands.
BOOST_AUTO_TEST_CASE(SurfaceExtentSpansLayersAndMaterial) {
  const Acts::MaterialSlab silicon =
      materialSlab(95.7f, 465.2f, 28.03f, 14.f, 8.28e-8f, 1.4f);

  BarrelEndcapDescription description;
  description.beamPipeRadius = 25.f;
  description.beamPipeMaterial = SurfaceMaterial{silicon};
  description.barrelRadii = {34.f, 99.f};
  description.barrelHalfLengthsZ = {244.f, 244.f};
  // the first barrel carries material out to where it ends, the second beyond
  // it, a service running past the last module
  description.barrelMaterials = {
      SurfaceMaterial{BandComposition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f},
                      {244.f},
                      {95.7f},
                      {465.2f}},
      SurfaceMaterial{BandComposition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f},
                      {244.f, 400.f},
                      {95.7f, 95.7f},
                      {465.2f, 465.2f}}};
  // a disc with a beam hole of nothing, one ring, and support past it
  description.discs = {
      {600.f,
       {{40.f, 190.f}},
       SurfaceMaterial{BandComposition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f},
                       {40.f, 190.f, 260.f},
                       {0.f, 95.7f, 95.7f},
                       {0.f, 465.2f, 465.2f}}}};
  description.passiveSurfaces = {
      {SurfaceShape::Cylinder, 60.f, 100.f, 300.f, SurfaceMaterial{silicon}}};

  const DetectorLayout layout = makeLayout(description);
  checkLayoutInvariants(layout);

  for (const DetectorSurface& surface : layout.surfaces) {
    if (surface.refCoord == 25.f) {
      // the beam pipe, which ends where the detector does
      BOOST_CHECK_EQUAL(surface.minBound, 0.f);
      BOOST_CHECK_EQUAL(surface.maxBound, 600.f);
    } else if (surface.refCoord == 60.f) {
      // a service, bounded by what it was declared with rather than by the
      // material it carries everywhere
      BOOST_CHECK_EQUAL(surface.minBound, 100.f);
      BOOST_CHECK_EQUAL(surface.maxBound, 300.f);
    } else if (surface.refCoord == 34.f) {
      // bounded by its modules, which its material stops with
      BOOST_CHECK_EQUAL(surface.minBound, 0.f);
      BOOST_CHECK_EQUAL(surface.maxBound, 244.f);
    } else if (surface.refCoord == 99.f) {
      // the service past the last module extends it
      BOOST_CHECK_EQUAL(surface.minBound, 0.f);
      BOOST_CHECK_EQUAL(surface.maxBound, 400.f);
    } else {
      // the disc: the ring bounds it below, the support above, and the hole
      // holds nothing to bound it with
      BOOST_CHECK_EQUAL(surface.minBound, 40.f);
      BOOST_CHECK_EQUAL(surface.maxBound, 260.f);
    }
  }

  // A surface of one material all over has no end, so a crossing of it is
  // always tested in full.
  BarrelEndcapDescription uniform;
  uniform.barrelRadii = {34.f};
  uniform.barrelHalfLengthsZ = {244.f};
  uniform.barrelMaterials = {SurfaceMaterial{silicon}};
  const DetectorLayout unbounded = makeLayout(uniform);
  BOOST_CHECK_EQUAL(unbounded.surfaces.front().maxBound,
                    std::numeric_limits<float>::infinity());
}

/// What every shipped description has to satisfy: the layout accounts for each
/// of its surfaces and rings, and the beam pipe is passive and in front.
BOOST_AUTO_TEST_CASE(ShippedLayouts) {
  const std::array descriptions{std::pair{genericDetectorPixelDescription(),
                                          makeGenericDetectorPixelLayout()}};

  for (const auto& [description, layout] : descriptions) {
    std::size_t rings = 0;
    for (const DiscDescription& disc : description.discs) {
      rings += disc.rings.size();
    }
    // beam pipe, the services, the barrel cylinders, and each disc on both
    // sides. A passive disc is mirrored, a passive cylinder is not.
    BOOST_CHECK_EQUAL(layout.surfaces.size(),
                      1u + passiveSurfaceCount(description) +
                          description.barrelRadii.size() +
                          2u * description.discs.size());
    BOOST_CHECK_EQUAL(layout.layers.size(), description.barrelRadii.size() *
                                                    description.barrelModules +
                                                2u * rings);
    checkLayoutInvariants(layout);

    BOOST_CHECK(layout.surfaces[0].layers.empty());
    BOOST_CHECK_LT(layout.surfaces[0].refCoord,
                   description.barrelRadii.front());
  }
}

/// A shipped description can be taken, modified and rebuilt: finer modules and
/// rings, and no beam pipe, which is legal and just leaves nothing in front of
/// the innermost layer. A modification that leaves it inconsistent throws
/// rather than defaulting the missing half of it.
BOOST_AUTO_TEST_CASE(DescriptionCanBeModifiedAndRebuilt) {
  BarrelEndcapDescription description = genericDetectorPixelDescription();
  description.barrelModules = 6;
  description.beamPipeRadius.reset();
  std::size_t rings = 0;
  for (DiscDescription& disc : description.discs) {
    disc.rings = subdivideRings(disc.rings, 2);
    rings += disc.rings.size();
  }
  const DetectorLayout layout = makeLayout(description);

  BOOST_CHECK_EQUAL(rings, 2u * 7u);
  BOOST_CHECK_EQUAL(layout.surfaces.size(),
                    passiveSurfaceCount(description) + 4u + 2u * 7u);
  BOOST_CHECK_EQUAL(layout.layers.size(), 4u * 6u + 2u * rings);
  // nothing sits where the beam pipe was; the services are further out
  for (const DetectorSurface& surface : layout.surfaces) {
    BOOST_CHECK_NE(surface.refCoord, 19.f);
  }
  checkLayoutInvariants(layout);

  description.barrelHalfLengthsZ.pop_back();
  BOOST_CHECK_THROW(makeLayout(description), std::invalid_argument);
}

/// `uniformRings` tiles a radial range exactly, and `subdivideRings` splits
/// each ring in place, so a gap between two survives.
BOOST_AUTO_TEST_CASE(RingHelpers) {
  const std::vector<RingBounds> uniform = uniformRings(30.f, 350.f, 4);
  BOOST_REQUIRE_EQUAL(uniform.size(), 4u);
  BOOST_CHECK_EQUAL(uniform.front().rMin, 30.f);
  BOOST_CHECK_EQUAL(uniform.back().rMax, 350.f);
  for (std::size_t m = 1; m < uniform.size(); ++m) {
    BOOST_CHECK_EQUAL(uniform[m].rMin, uniform[m - 1].rMax);
  }
  BOOST_CHECK_THROW(uniformRings(350.f, 30.f, 1), std::invalid_argument);
  BOOST_CHECK_THROW(uniformRings(30.f, 350.f, 0), std::invalid_argument);

  const std::array coarse{RingBounds{30.f, 50.f}, RingBounds{100.f, 140.f}};
  const std::vector<RingBounds> split = subdivideRings(coarse, 2);
  BOOST_REQUIRE_EQUAL(split.size(), 4u);
  BOOST_CHECK_EQUAL(split[0].rMin, 30.f);
  BOOST_CHECK_EQUAL(split[1].rMax, 50.f);
  BOOST_CHECK_EQUAL(split[2].rMin, 100.f);
  BOOST_CHECK_EQUAL(split[3].rMax, 140.f);
  // the gap between the two rings is still there
  BOOST_CHECK_LT(split[1].rMax, split[2].rMin);
  // and asking for one part leaves them alone
  BOOST_CHECK_EQUAL(subdivideRings(coarse, 1).size(), 2u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
