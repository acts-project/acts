// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Tests/CommonHelpers/DataDirectory.hpp"

#include <fstream>

class DummyDecorator : public Acts::IVolumeMaterialJsonDecorator {
 public:
  void decorate([[maybe_unused]] const Acts::ISurfaceMaterial &material,
                [[maybe_unused]] nlohmann::json &json) const override{};

  void decorate([[maybe_unused]] const Acts::IVolumeMaterial &material,
                [[maybe_unused]] nlohmann::json &json) const override{};
};

BOOST_AUTO_TEST_SUITE(MaterialMapJsonConverter)

BOOST_AUTO_TEST_CASE(RoundtripFromFile) {
  // read reference map from file
  std::ifstream refFile(Acts::Test::getDataPath("material-map.json"));
  nlohmann::json refJson;
  refFile >> refJson;

  DummyDecorator decorator;
  // convert to the material map and back again
  Acts::MaterialMapJsonConverter::Config converterCfg;
  Acts::MaterialMapJsonConverter converter(converterCfg, Acts::Logging::INFO);
  auto materialMap = converter.jsonToMaterialMaps(refJson);
  nlohmann::json encodedJson =
      converter.materialMapsToJson(materialMap, &decorator);

  // verify identical encoded JSON values
  BOOST_CHECK_EQUAL(refJson, encodedJson);
}

BOOST_AUTO_TEST_SUITE_END()
