// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsFatras/Json/DataDirectory.hpp"
#include "ActsFatras/Json/DetectorDescriptionJsonConverter.hpp"
#include "ActsFatras/Json/EventConfigJsonConverter.hpp"
#include "ActsFatras/Synthetic/DetectorLayout.hpp"
#include "ActsFatras/Synthetic/EventConfig.hpp"
#include "ActsFatras/Synthetic/EventGenerator.hpp"
#include "ActsFatras/Synthetic/SyntheticEvent.hpp"

#include <cstddef>
#include <string>
#include <vector>

using namespace ActsFatras::Synthetic;

namespace ActsTests {

namespace {

/// One of the detectors that ships in `Fatras/data`, read back the way a caller
/// reads it: the description, then the material onto it.
/// @param name the detector, e.g. `itk`
/// @return its description, decorated
DetectorDescription readShipped(const std::string& name) {
  DetectorDescription description =
      readDetectorDescription(ActsFatras::dataPath(name + "-description.json"));
  decorate(description, readMaterialDecoration(
                            ActsFatras::dataPath(name + "-material.json")));
  return description;
}

/// @param subsystem the subsystem to count
/// @return how many cylinders its barrels hold between them
std::size_t cylinderCount(const SubsystemDescription& subsystem) {
  std::size_t count = 0;
  for (const BarrelDescription& barrel : subsystem.barrels) {
    count += barrel.cylinders.size();
  }
  return count;
}

/// @param subsystem the subsystem to count
/// @return how many discs its endcaps hold on one side, and how many rings
std::pair<std::size_t, std::size_t> discAndRingCount(
    const SubsystemDescription& subsystem) {
  std::size_t discs = 0;
  std::size_t rings = 0;
  for (const EndcapDescription& endcap : subsystem.endcaps) {
    discs += endcap.discs.size();
    for (const DiscDescription& disc : endcap.discs) {
      rings += disc.rings.size();
    }
  }
  return {discs, rings};
}

/// Every layer of a detector that is meant to hold material holds some, and
/// nothing holds a negative amount of it.
/// @param layout the detector to check
/// @return how many of its surfaces carry material
std::size_t checkMaterial(const DetectorLayout& layout) {
  std::size_t carrying = 0;
  for (const DetectorSurface& surface : layout.surfaces) {
    if (surface.material.bands.empty()) {
      continue;
    }
    ++carrying;
    bool anything = false;
    for (const Acts::MaterialSlab& band : surface.material.bands) {
      BOOST_CHECK_GE(band.thicknessInX0(), 0.f);
      BOOST_CHECK_GE(band.thicknessInL0(), 0.f);
      anything = anything || !band.isVacuum();
    }
    BOOST_CHECK(anything);
  }
  return carrying;
}

/// A detector and a configuration make an event a track can be found in.
/// @param layout the detector to generate on
/// @param config the configuration to generate with
void checkGeneratesAnEvent(const DetectorLayout& layout,
                           const EventConfig& config) {
  const Event event = generateEvent(layout, config);
  BOOST_REQUIRE(!event.particles.empty());
  BOOST_REQUIRE(!event.spacePoints.empty());

  // A primary of this detector leaves of the order of ten space points. Loose
  // on purpose: the tuning is checked against a full simulation elsewhere, and
  // what is being checked here is that the files describe a detector a track
  // can cross at all.
  std::size_t primaries = 0;
  std::size_t primaryHits = 0;
  for (const GeneratedParticle& particle : event.particles) {
    if (!particle.primary()) {
      continue;
    }
    ++primaries;
    primaryHits += particle.numHits;
  }
  BOOST_REQUIRE_GT(primaries, 0u);
  const double hitsPerPrimary =
      static_cast<double>(primaryHits) / static_cast<double>(primaries);
  BOOST_CHECK_GT(hitsPerPrimary, 2.);
  BOOST_CHECK_LT(hitsPerPrimary, 20.);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticShippedDetectorsSuite)

/// The ITk pixel description as it ships: the ATLAS layout it was transcribed
/// from, down to the number of rings, and the beam pipe outside the pixels
/// because that is where the real geometry keeps it.
BOOST_AUTO_TEST_CASE(ItkPixel) {
  const DetectorDescription description = readShipped("itk");

  BOOST_CHECK_EQUAL(description.escapeRadius, 1000.f);
  BOOST_CHECK_EQUAL(description.escapeHalfZ, 3050.f);

  // the beam pipe belongs to the detector rather than to the pixels
  BOOST_REQUIRE_EQUAL(description.passives.size(), 1u);
  BOOST_CHECK_EQUAL(description.passives.front().refCoord, 25.f);
  BOOST_CHECK(description.passives.front().shape == SurfaceShape::Cylinder);

  BOOST_REQUIRE_EQUAL(description.subsystems.size(), 1u);
  const SubsystemDescription& pixels = description.subsystems.front();
  BOOST_CHECK_EQUAL(pixels.name, "itk-pixel");
  // five stave layers, and seventy-five discs per side carrying ninety-five
  // rings: the endcap is the staggered thing it is, not a stack of uniform
  // discs
  BOOST_CHECK_EQUAL(cylinderCount(pixels), 5u);
  const auto [discs, rings] = discAndRingCount(pixels);
  BOOST_CHECK_EQUAL(discs, 75u);
  BOOST_CHECK_EQUAL(rings, 95u);
  // the services, whose extent matters as much as their amount
  BOOST_CHECK_EQUAL(pixels.passives.size(), 4u);

  const DetectorLayout layout = makeLayout(description);
  // beam pipe, three service cylinders and a service disc on both sides, the
  // staves, and every disc on both sides
  BOOST_CHECK_EQUAL(layout.surfaces.size(), 1u + 5u + 5u + 2u * 75u);
  BOOST_CHECK_EQUAL(layout.layers.size(), 5u + 2u * 95u);
  BOOST_REQUIRE_EQUAL(layout.subsystems.size(), 1u);
  BOOST_CHECK_EQUAL(layout.subsystems.front(), "itk-pixel");
  // every sensitive surface and every passive carries material
  BOOST_CHECK_EQUAL(checkMaterial(layout), layout.surfaces.size());
}

/// The configuration fitted to it, which has to be read as a whole.
BOOST_AUTO_TEST_CASE(ItkPixelConfiguration) {
  const EventConfig config =
      readEventConfig(ActsFatras::dataPath("itk-ttbar-pu200.json"));

  BOOST_CHECK_EQUAL(config.generation.pileup, 200u);
  BOOST_CHECK_CLOSE(config.generation.chargedPerUnitRapidity, 6.15f, 1e-3);
  // only a curler gets past two and a half turns, which the reference's soft
  // primaries do
  BOOST_CHECK_CLOSE(config.simulation.propagation.maxTurns, 3.f, 1e-3);
  BOOST_CHECK_CLOSE(config.simulation.secondaries.electronRate, 3.202f, 1e-3);
  BOOST_CHECK_CLOSE(config.simulation.secondaries.nuclearRate, 7.483f, 1e-3);
  BOOST_CHECK(config.simulation.material.energyLossModel ==
              EnergyLossModel::Mode);
}

/// The Open Data Detector pixels, read off the geometry ACTS builds for it and
/// kept as data so that the configuration fitted to them works without DD4hep.
BOOST_AUTO_TEST_CASE(OpenDataDetectorPixel) {
  const DetectorDescription description = readShipped("odd");

  BOOST_CHECK_EQUAL(description.escapeRadius, 1100.f);
  BOOST_CHECK_EQUAL(description.escapeHalfZ, 3000.f);

  // the beryllium beam pipe, which sits in a volume of its own that a
  // pixel-only selector excludes
  BOOST_REQUIRE_EQUAL(description.passives.size(), 1u);
  BOOST_CHECK_EQUAL(description.passives.front().refCoord, 24.f);

  BOOST_REQUIRE_EQUAL(description.subsystems.size(), 1u);
  const SubsystemDescription& pixels = description.subsystems.front();
  BOOST_CHECK_EQUAL(pixels.name, "odd-pixel");
  BOOST_CHECK_EQUAL(cylinderCount(pixels), 4u);
  // Seven disc positions a side, each read off the geometry as the two z planes
  // seven millimetres apart that its modules alternate between, so fourteen
  // discs of one ring each. That is what a track crossing them meets.
  const auto [discs, rings] = discAndRingCount(pixels);
  BOOST_CHECK_EQUAL(discs, 14u);
  BOOST_CHECK_EQUAL(rings, 14u);
  const std::vector<DiscDescription>& oddDiscs = pixels.endcaps.front().discs;
  BOOST_CHECK_CLOSE(oddDiscs[1].absZ - oddDiscs[0].absZ, 7.f, 1e-2);
  // the two service discs the reduction found away from any layer, the outer
  // one beyond the last disc entirely
  BOOST_REQUIRE_EQUAL(pixels.passives.size(), 2u);
  BOOST_CHECK_CLOSE(pixels.passives[0].refCoord, 559.f, 1e-3);
  BOOST_CHECK_CLOSE(pixels.passives[1].refCoord, 1950.f, 1e-3);
  BOOST_CHECK_GT(pixels.passives[1].refCoord, oddDiscs.back().absZ);

  const DetectorLayout layout = makeLayout(description);
  // beam pipe, the two service discs on both sides, the barrels, and the discs
  // on both sides
  BOOST_CHECK_EQUAL(layout.surfaces.size(), 1u + 2u * 2u + 4u + 2u * 14u);
  BOOST_CHECK_EQUAL(layout.layers.size(), 4u + 2u * 14u);
  BOOST_CHECK_EQUAL(checkMaterial(layout), layout.surfaces.size());
}

/// The configuration fitted to it, against ColliderML rather than an ATLAS
/// dump.
BOOST_AUTO_TEST_CASE(OpenDataDetectorConfiguration) {
  const EventConfig config =
      readEventConfig(ActsFatras::dataPath("odd-ttbar-pu200.json"));

  BOOST_CHECK_EQUAL(config.generation.pileup, 200u);
  BOOST_CHECK_GT(config.simulation.secondaries.electronRate, 0.f);
  BOOST_CHECK_GT(config.simulation.secondaries.nuclearRate, 0.f);
  // its own detector's resolution, not the ITk's
  BOOST_CHECK_GT(config.simulation.measurement.positionSmearing, 0.f);
}

/// A detector and the configuration fitted to it generate an event: what the
/// two files are for, and the one check that exercises every number in them.
BOOST_AUTO_TEST_CASE(ShippedDetectorsGenerateAnEvent) {
  for (const std::string& name : {std::string{"itk"}, std::string{"odd"}}) {
    const DetectorLayout layout = makeLayout(readShipped(name));
    EventConfig config =
        readEventConfig(ActsFatras::dataPath(name + "-ttbar-pu200.json"));
    // a few interactions rather than the two hundred it is fitted at, this
    // being a test of the files and not of the tuning
    config.generation.pileup = 2;
    checkGeneratesAnEvent(layout, config);
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
