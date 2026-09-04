// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Json/IVolumeMaterialJsonDecorator.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
#include "ActsTests/CommonHelpers/DataDirectory.hpp"

#include <fstream>
#include <memory>

#include <nlohmann/json.hpp>

namespace Acts {
class IVolumeMaterial;
}  // namespace Acts

using namespace Acts;

class DummyDecorator : public IVolumeMaterialJsonDecorator {
 public:
  void decorate([[maybe_unused]] const ISurfaceMaterial &material,
                [[maybe_unused]] nlohmann::json &json) const override {};

  void decorate([[maybe_unused]] const IVolumeMaterial &material,
                [[maybe_unused]] nlohmann::json &json) const override {};
};

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(JsonSuite)

BOOST_AUTO_TEST_CASE(RoundtripFromFile) {
  // read reference map from file
  std::ifstream refFile(ActsTests::getDataPath("material-map.json"));
  nlohmann::json refJson;
  refFile >> refJson;

  DummyDecorator decorator;
  // convert to the material map and back again
  MaterialMapJsonConverter::Config converterCfg;
  MaterialMapJsonConverter converter(converterCfg, Logging::INFO);
  auto materialMap = converter.jsonToMaterialMaps(refJson);
  nlohmann::json encodedJson =
      converter.materialMapsToJson(materialMap, &decorator);

  // verify identical encoded JSON values
  BOOST_CHECK_EQUAL(refJson, encodedJson);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
