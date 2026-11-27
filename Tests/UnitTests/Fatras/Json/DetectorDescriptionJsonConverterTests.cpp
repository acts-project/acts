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
#include "ActsFatras/Synthetic/DetectorLayout.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using namespace ActsFatras::Synthetic;

namespace ActsTests {

namespace {

/// A description numbered the way `decorate` and `makeLayout` number it, which
/// is what a description that has been through a file comes back as.
/// @param description the detector to number
/// @return a numbered copy of it
DetectorDescription numbered(DetectorDescription description) {
  assignLayerIndices(description);
  return description;
}

/// Two layouts agree where it matters: the same surfaces in the same order,
/// each bounded the same way and made of the same thing.
/// @param actual the layout to check
/// @param expected what it has to come to
void checkSameLayout(const DetectorLayout& actual,
                     const DetectorLayout& expected) {
  BOOST_REQUIRE_EQUAL(actual.surfaces.size(), expected.surfaces.size());
  BOOST_REQUIRE_EQUAL(actual.layers.size(), expected.layers.size());
  BOOST_CHECK(actual.subsystems == expected.subsystems);
  BOOST_CHECK_EQUAL(actual.escapeRadius, expected.escapeRadius);
  BOOST_CHECK_EQUAL(actual.escapeHalfZ, expected.escapeHalfZ);

  for (std::size_t s = 0; s < expected.surfaces.size(); ++s) {
    const DetectorSurface& was = expected.surfaces[s];
    const DetectorSurface& is = actual.surfaces[s];
    BOOST_CHECK_EQUAL(is.refCoord, was.refCoord);
    BOOST_CHECK_EQUAL(is.minBound, was.minBound);
    BOOST_CHECK_EQUAL(is.maxBound, was.maxBound);
    BOOST_CHECK_EQUAL(is.overlapProbability, was.overlapProbability);
    BOOST_CHECK_EQUAL(is.overlapOffset, was.overlapOffset);
    BOOST_CHECK(is.shape == was.shape);
    BOOST_CHECK(is.layers == was.layers);

    BOOST_REQUIRE_EQUAL(is.material.bands.size(), was.material.bands.size());
    BOOST_CHECK(is.material.bounds == was.material.bounds);
    for (std::size_t k = 0; k < was.material.bands.size(); ++k) {
      BOOST_CHECK_EQUAL(is.material.bands[k].thickness(),
                        was.material.bands[k].thickness());
      BOOST_CHECK_EQUAL(is.material.bands[k].thicknessInX0(),
                        was.material.bands[k].thicknessInX0());
      BOOST_CHECK_EQUAL(is.material.bands[k].thicknessInL0(),
                        was.material.bands[k].thicknessInL0());
    }
  }
  for (std::size_t l = 0; l < expected.layers.size(); ++l) {
    BOOST_CHECK_EQUAL(actual.layers[l].refCoord, expected.layers[l].refCoord);
    BOOST_CHECK_EQUAL(actual.layers[l].minBound, expected.layers[l].minBound);
    BOOST_CHECK_EQUAL(actual.layers[l].maxBound, expected.layers[l].maxBound);
    BOOST_CHECK_EQUAL(actual.layers[l].layer, expected.layers[l].layer);
    BOOST_CHECK_EQUAL(actual.layers[l].moduleIndex,
                      expected.layers[l].moduleIndex);
    BOOST_CHECK_EQUAL(actual.layers[l].subsystem, expected.layers[l].subsystem);
  }
}

}  // namespace

BOOST_AUTO_TEST_SUITE(SyntheticDescriptionJsonSuite)

/// A description survives a round trip through JSON: the same file comes back
/// out, and it expands to the same detector.
BOOST_AUTO_TEST_CASE(DescriptionRoundTrip) {
  const DetectorDescription original = genericDetectorPixelDescription();

  const nlohmann::json written = original;
  const DetectorDescription read = written.get<DetectorDescription>();
  const nlohmann::json rewritten = read;
  BOOST_CHECK(written == rewritten);

  checkSameLayout(makeLayout(read), makeLayout(original));
}

/// The same through a file, header and all. A file that is not one of these is
/// refused rather than read as an empty detector.
BOOST_AUTO_TEST_CASE(DescriptionRoundTripsThroughAFile) {
  const DetectorDescription original = genericDetectorPixelDescription();
  const std::filesystem::path path = std::filesystem::temp_directory_path() /
                                     "acts-synthetic-description.json";

  writeDetectorDescription(path, original);
  const DetectorDescription read = readDetectorDescription(path);
  BOOST_CHECK(nlohmann::json(read) == nlohmann::json(original));
  checkSameLayout(makeLayout(read), makeLayout(original));

  // the header says what kind of file it is, and it is checked
  const std::filesystem::path other =
      std::filesystem::temp_directory_path() / "acts-synthetic-not-one.json";
  {
    std::ofstream file(other);
    file << nlohmann::json(original).dump(2) << std::endl;
  }
  BOOST_CHECK_THROW(readDetectorDescription(other), std::runtime_error);
  BOOST_CHECK_THROW(
      readDetectorDescription(std::filesystem::temp_directory_path() /
                              "acts-synthetic-does-not-exist.json"),
      std::runtime_error);

  std::filesystem::remove(path);
  std::filesystem::remove(other);
}

/// Where a detector is and what it is made of are two files: taking the
/// material off a description and putting it back gives the description back.
BOOST_AUTO_TEST_CASE(MaterialSplitsOffAndGoesBackOn) {
  const DetectorDescription original = genericDetectorPixelDescription();
  const MaterialDecoration decoration = extractMaterial(original);
  BOOST_CHECK(!decoration.empty());

  DetectorDescription bare = original;
  stripMaterial(bare);
  // a description of where the layers are and nothing about what they hold
  BOOST_CHECK(nlohmann::json(bare).dump().find("\"material\"") ==
              std::string::npos);
  BOOST_CHECK(extractMaterial(bare).empty());
  // and the geometry is untouched by the stripping
  const DetectorLayout stripped = makeLayout(bare);
  BOOST_CHECK_EQUAL(stripped.surfaces.size(),
                    makeLayout(original).surfaces.size());
  for (const DetectorSurface& surface : stripped.surfaces) {
    BOOST_CHECK(surface.material.bands.empty());
  }

  // the decoration goes through a file of its own
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "acts-synthetic-material.json";
  writeMaterialDecoration(path, decoration);
  const MaterialDecoration read = readMaterialDecoration(path);
  BOOST_CHECK(nlohmann::json(read) == nlohmann::json(decoration));
  std::filesystem::remove(path);

  decorate(bare, read);
  BOOST_CHECK(nlohmann::json(bare) == nlohmann::json(numbered(original)));
  checkSameLayout(makeLayout(bare), makeLayout(original));
}

/// A key that names no layer is a material file left behind by a description
/// that has been renumbered, which is the one thing the explicit layer indices
/// exist to catch.
BOOST_AUTO_TEST_CASE(DecorationRejectsAStaleKey) {
  DetectorDescription description = genericDetectorPixelDescription();
  const MaterialDecoration decoration = extractMaterial(description);
  BOOST_REQUIRE(!decoration.empty());

  MaterialDecoration stale = decoration;
  stale.front().layer.layer = 99;
  BOOST_CHECK_THROW(decorate(description, stale), std::invalid_argument);

  MaterialDecoration renamed = decoration;
  renamed.front().layer.subsystem = "generic-strip";
  BOOST_CHECK_THROW(decorate(description, renamed), std::invalid_argument);

  // and the message says which layer it could not find
  try {
    decorate(description, stale);
  } catch (const std::invalid_argument& error) {
    BOOST_CHECK_NE(std::string(error.what()).find("layer 99"),
                   std::string::npos);
  }
}

/// An unbounded surface has no number for where it ends, so it says so.
BOOST_AUTO_TEST_CASE(InfiniteBoundsSurviveTheRoundTrip) {
  DetectorDescription description;
  description.subsystems = {
      SubsystemDescription{.name = "pixel",
                           .passives = {PassiveSurfaceDescription{
                               .shape = SurfaceShape::Cylinder,
                               .refCoord = 124.f,
                               .material = SurfaceMaterial{Acts::MaterialSlab{
                                   Acts::Material::fromMolarDensity(
                                       95.7f, 465.2f, 28.f, 14.f, 8.28e-8f),
                                   1.4f}}}}}};

  const nlohmann::json written = description;
  const std::string dumped = written.dump();
  // both the surface's own end and the outer edge of its one band
  BOOST_CHECK_NE(dumped.find("\"inf\""), std::string::npos);

  const DetectorDescription read = written.get<DetectorDescription>();
  const PassiveSurfaceDescription& passive =
      read.subsystems.front().passives.front();
  BOOST_CHECK_EQUAL(passive.maxBound, std::numeric_limits<float>::infinity());
  BOOST_CHECK_EQUAL(passive.material.bounds.back(),
                    std::numeric_limits<float>::infinity());
  checkSameLayout(makeLayout(read), makeLayout(description));
}

/// A number that every band agrees on is written once. Which is not a saving
/// worth having on its own -- it is what keeps a surface banded only in how
/// much material it holds readable, its composition being stated where it
/// belongs.
BOOST_AUTO_TEST_CASE(BandNumbersCollapseWhereTheyAgree) {
  const BandComposition composition{28.03f, 14.f, 8.28e-8f * 95.7f, 1.4f};
  const SurfaceMaterial oneComposition{
      composition, {30.f, 50.f, 80.f}, {95.7f, 957.f}, {465.2f, 1914.f}};
  const nlohmann::json single = oneComposition;
  // the two lengths differ band by band, so they are spelled out
  BOOST_CHECK(single.at("x0").is_array());
  BOOST_CHECK(single.at("l0").is_array());
  // what it is made of does not, so it is stated once -- and the one number
  // that stands for the density is the product with the radiation length, the
  // density following that length band by band
  BOOST_CHECK(single.at("ar").is_number());
  BOOST_CHECK(single.at("z").is_number());
  BOOST_CHECK(single.at("molarDensityX0").is_number());
  BOOST_CHECK(single.at("thickness").is_number());
  BOOST_CHECK_EQUAL(single.at("ar").get<float>(), 28.03f);
  BOOST_CHECK_CLOSE(single.at("molarDensityX0").get<float>(), 8.28e-8f * 95.7f,
                    1e-3);

  // A band of nothing has no composition to disagree with, so it does not force
  // the long form on the bands that do hold something. Its radiation length
  // does not collapse, being how the file says the band holds nothing at all.
  const SurfaceMaterial withAGap{composition,
                                 {30.f, 50.f, 80.f, 120.f},
                                 {95.7f, 0.f, 95.7f},
                                 {465.2f, 0.f, 465.2f}};
  const nlohmann::json gapped = withAGap;
  BOOST_CHECK(gapped.at("ar").is_number());
  BOOST_CHECK(gapped.at("z").is_number());
  BOOST_CHECK(gapped.at("molarDensityX0").is_number());
  BOOST_REQUIRE(gapped.at("x0").is_array());
  BOOST_CHECK_EQUAL(gapped.at("x0").at(1).get<float>(), 0.f);
  BOOST_CHECK(withAGap.bands[1].isVacuum());

  for (const SurfaceMaterial& material : {oneComposition, withAGap}) {
    const auto read = nlohmann::json(material).get<SurfaceMaterial>();
    BOOST_REQUIRE_EQUAL(read.bands.size(), material.bands.size());
    BOOST_CHECK(read.bounds == material.bounds);
    for (std::size_t k = 0; k < material.bands.size(); ++k) {
      BOOST_CHECK_EQUAL(read.bands[k].thicknessInX0(),
                        material.bands[k].thicknessInX0());
      BOOST_CHECK_EQUAL(read.bands[k].thicknessInL0(),
                        material.bands[k].thicknessInL0());
      BOOST_CHECK_EQUAL(read.bands[k].thickness(),
                        material.bands[k].thickness());
    }
  }

  // a per-band field of the wrong length has lost a number, and taking it would
  // shift every band along the surface
  nlohmann::json broken = single;
  broken["x0"] = std::vector<float>{95.7f};
  BOOST_CHECK_THROW(broken.get<SurfaceMaterial>(), std::invalid_argument);
}

/// The shipped Generic files say the same thing as the C++ description they
/// were written from, which is what pins the format against a detector that
/// exists in both forms.
BOOST_AUTO_TEST_CASE(ShippedGenericFilesMatchTheCode) {
  DetectorDescription fromFile = readDetectorDescription(
      ActsFatras::dataPath("generic-pixel-description.json"));
  decorate(fromFile, readMaterialDecoration(
                         ActsFatras::dataPath("generic-pixel-material.json")));

  const DetectorDescription fromCode =
      numbered(genericDetectorPixelDescription());
  BOOST_CHECK(nlohmann::json(fromFile) == nlohmann::json(fromCode));
  checkSameLayout(makeLayout(fromFile), makeLayout(fromCode));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
