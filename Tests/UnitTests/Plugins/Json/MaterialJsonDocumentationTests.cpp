// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsTests/CommonHelpers/PredefinedMaterials.hpp"

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

using namespace Acts;
using namespace Acts::UnitLiterals;

namespace ActsTests {

namespace {

/// The one real material of the example, reused wherever material is needed
MaterialSlab exampleSlab() {
  return {makeSilicon(), 0.15_mm};
}

/// The store the index grid points into. Vacuum sits at index 0 so the example
/// shows both an empty slab and an index that is not the bin number.
std::vector<MaterialSlab> exampleSlabStore() {
  return {MaterialSlab::Nothing(), exampleSlab()};
}

/// The material of the documented example.
///
/// This is not a coverage test -- that the encoder handles every material type
/// is checked by SurfaceMaterialJsonConverterTests. The example only has to
/// show the structure, so it carries the three surface payloads that differ in
/// shape (a plain slab, a binned matrix, a grid) plus one volume entry, and
/// every count is the smallest the format allows. The result is included
/// verbatim in the documentation, so every line of it is a line someone reads.
TrackingGeometryMaterial exampleMaterialMaps() {
  SurfaceMaterialMaps surfaces;
  VolumeMaterialMaps volumes;

  const GeometryIdentifier volume1 = GeometryIdentifier().withVolume(1);
  const GeometryIdentifier volume2 = GeometryIdentifier().withVolume(2);

  // A single slab covering the whole surface
  surfaces[volume1.withBoundary(1)] =
      std::make_shared<const HomogeneousSurfaceMaterial>(
          exampleSlab(), 1., MappingType::PreMapping);

  // The classic binned material: a BinUtility plus a slab matrix. The second
  // bin is left empty to show how an uncovered bin is written.
  BinUtility binUtility(2, -100., 100., open, AxisDirection::AxisZ);
  MaterialSlabMatrix matrix{
      MaterialSlabVector{exampleSlab(), MaterialSlab::Nothing()}};
  surfaces[volume1.withLayer(2).withApproach(1)] =
      std::make_shared<const BinnedSurfaceMaterial>(binUtility,
                                                    std::move(matrix));

  // Grid material, bins index into a slab store. The other two storage
  // backends, and the proto materials, differ from this only in the keys
  // listed in the docs. The grid is always two dimensional.
  {
    auto axis0 =
        IAxis::createEquidistant(AxisBoundaryType::Bound, -100., 100., 1);
    auto axis1 =
        IAxis::createEquidistant(AxisBoundaryType::Bound, -100., 100., 1);
    surfaces[volume2.withLayer(4)] = GridSurfaceMaterial::createIndexed(
        *axis0, *axis1, exampleSlabStore(),
        std::vector<std::vector<std::size_t>>{std::vector<std::size_t>{1u}});
  }

  volumes[volume1] =
      std::make_shared<const HomogeneousVolumeMaterial>(makeSilicon());

  return {std::move(surfaces), std::move(volumes)};
}

/// The documented example file to check against, passed on the command line by
/// ctest so the test does not have to know where the source tree lives.
std::string exampleFilePath() {
  const auto& master = boost::unit_test::framework::master_test_suite();
  BOOST_REQUIRE_MESSAGE(master.argc == 2,
                        "Expected exactly one argument, the example file path");
  return master.argv[1];
}

}  // namespace

BOOST_AUTO_TEST_SUITE(JsonSuite)

/// Writes the material map example that the JSON plugin documentation shows,
/// and fails if the committed file no longer matches what the converter
/// produces. Set ACTS_UPDATE_DOC_EXAMPLES=1 to refresh the file in place.
BOOST_AUTO_TEST_CASE(MaterialMapDocumentationExample) {
  MaterialMapJsonConverter::Config converterCfg;
  MaterialMapJsonConverter converter(converterCfg, Logging::WARNING);

  nlohmann::json jMap = converter.materialMapsToJson(exampleMaterialMaps());

  // The documented example has to be something the plugin can read back
  nlohmann::json jRoundTrip =
      converter.materialMapsToJson(converter.jsonToMaterialMaps(jMap));
  BOOST_CHECK_EQUAL(jMap, jRoundTrip);

  const std::string path = exampleFilePath();
  // Four spaces and a trailing newline, so the file survives the
  // repository's json formatting hook unchanged
  const std::string encoded = jMap.dump(4) + "\n";

  if (const char* update = std::getenv("ACTS_UPDATE_DOC_EXAMPLES");
      update != nullptr && std::string(update) != "0") {
    std::ofstream out(path);
    BOOST_REQUIRE_MESSAGE(out.is_open(), "Cannot write " << path);
    out << encoded;
    return;
  }

  std::ifstream in(path);
  BOOST_REQUIRE_MESSAGE(in.is_open(), "Cannot read " << path);
  const std::string committed{std::istreambuf_iterator<char>(in),
                              std::istreambuf_iterator<char>()};

  BOOST_CHECK_MESSAGE(
      committed == encoded,
      path << " is out of date, regenerate it by running this test with "
              "ACTS_UPDATE_DOC_EXAMPLES=1");
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
