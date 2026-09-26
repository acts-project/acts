// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

using namespace Acts;
using Converter = TrackingGeometryMaterialJsonConverter;

namespace ActsTests {
namespace {
std::filesystem::path examples() {
  const auto& master = boost::unit_test::framework::master_test_suite();
  BOOST_REQUIRE_EQUAL(master.argc, 2);
  return master.argv[1];
}
nlohmann::json fixture(const char* name) {
  return detail::readJsonFile(examples() / name);
}

class CustomMaterial final : public ISurfaceMaterial {
 public:
  using ISurfaceMaterial::materialSlab;
  explicit CustomMaterial(double thickness)
      : m_slab(MaterialSlab::Vacuum(static_cast<float>(thickness))) {}
  CustomMaterial& scale(double factor) override {
    m_slab.scaleThickness(static_cast<float>(factor));
    return *this;
  }
  const MaterialSlab& materialSlab(
      const Vector2& /*localPosition*/) const override {
    return m_slab;
  }
  std::vector<AxisDirection> localAxisDirections() const override { return {}; }
  std::ostream& toStream(std::ostream& out) const override {
    return out << "custom";
  }

 private:
  MaterialSlab m_slab;
};
nlohmann::json encodeCustom(const CustomMaterial& m,
                            Converter::EncodeContext& /*context*/) {
  return {{"kind", "application-custom"},
          {"thickness", m.materialSlab(Vector2::Zero()).thickness()}};
}
std::unique_ptr<const ISurfaceMaterial> decodeCustom(
    const nlohmann::json& j, const Converter::DecodeContext& /*context*/) {
  return std::make_unique<CustomMaterial>(j.at("thickness").get<double>());
}
}  // namespace

BOOST_AUTO_TEST_SUITE(JsonSuite)
BOOST_AUTO_TEST_CASE(MaterialDocumentExamples) {
  Converter converter;
  for (const auto* name : {"minimal.json", "surfaces.json", "templates.json"}) {
    BOOST_TEST_CONTEXT(name) {
      const auto decoded = converter.fromFile(examples() / name);
      const auto encoded = converter.toJson(decoded);
      BOOST_CHECK(encoded == converter.toJson(converter.fromJson(encoded)));
      // Optional export for offline schema validation of actual codec output.
      if (const char* output =
              std::getenv("ACTS_MATERIAL_DOCUMENT_OUTPUT_DIR")) {
        std::filesystem::create_directories(output);
        std::ofstream(std::filesystem::path(output) / name)
            << encoded.dump(4) << '\n';
      }
    }
  }
  const auto surfaces = converter.fromJson(fixture("surfaces.json"));
  BOOST_REQUIRE(surfaces.description());
  BOOST_CHECK_EQUAL(*surfaces.description(), fixture("surfaces.json")
                                                 .at("header")
                                                 .at("description")
                                                 .get<std::string>());
  const auto* a = dynamic_cast<const GridSurfaceMaterial*>(
      surfaces.keyedSurfaces.at("barrel/support-a").material.get());
  const auto* b = dynamic_cast<const GridSurfaceMaterial*>(
      surfaces.keyedSurfaces.at("barrel/support-b").material.get());
  BOOST_REQUIRE(a != nullptr);
  BOOST_REQUIRE(b != nullptr);
  BOOST_CHECK(
      std::get<GridSurfaceMaterial::GloballyIndexed>(a->storage()).material ==
      std::get<GridSurfaceMaterial::GloballyIndexed>(b->storage()).material);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentDescription) {
  Converter converter;
  TrackingGeometryMaterial material;
  BOOST_CHECK(!material.description());
  BOOST_CHECK(!converter.toJson(material).at("header").contains("description"));
  material.setDescription("Mapped detector material");
  auto copy = material;
  copy.setDescription("");
  BOOST_CHECK_EQUAL(*material.description(), "Mapped detector material");
  BOOST_REQUIRE(copy.description());
  BOOST_CHECK(copy.description()->empty());
  const auto encoded = converter.toJson(material);
  BOOST_CHECK_EQUAL(*converter.fromJson(encoded).description(),
                    "Mapped detector material");
  BOOST_REQUIRE(converter.fromJson(converter.toJson(copy)).description());
  copy.setDescription(std::nullopt);
  BOOST_CHECK(!converter.toJson(copy).at("header").contains("description"));
  for (const nlohmann::json& invalid :
       {nlohmann::json(nullptr), nlohmann::json(42), nlohmann::json::object(),
        nlohmann::json::array()}) {
    auto malformed = encoded;
    malformed["header"]["description"] = invalid;
    BOOST_CHECK_THROW(converter.fromJson(malformed), nlohmann::json::exception);
  }
}

BOOST_AUTO_TEST_CASE(MaterialDocumentPhysicalPropertiesAndStorage) {
  Converter converter;
  const auto silicon = Material::fromMolarDensity(90.f, 450.f, 28.f, 14.f,
                                                  0.001f, 0.02f, 0.0003f);
  const MaterialSlab slab(silicon, 0.25f);
  TrackingGeometryMaterial source;
  const auto id = GeometryIdentifier()
                      .withVolume(255)
                      .withBoundary(255)
                      .withLayer(4095)
                      .withApproach(255)
                      .withSensitive(1048575)
                      .withExtra(255);
  source.surfaceMaterials[id] = std::make_shared<HomogeneousSurfaceMaterial>(
      slab, 0.25, MappingType::Sensor);
  source.setDescription("roundtrip");
  const auto recovered = converter.fromJson(converter.toJson(source));
  BOOST_CHECK(source.description() == recovered.description());
  const auto& material = *recovered.surfaceMaterials.at(id);
  BOOST_CHECK(material.materialSlab(Vector2::Zero()) == slab);
  BOOST_CHECK_EQUAL(
      material.factor(Direction::Negative(), MaterialUpdateMode::PreUpdate),
      0.25);
  BOOST_CHECK(material.mappingType() == MappingType::Sensor);

  // Unequal dimensions and unique values expose accidental native-order copies,
  // including all guard cells rather than only regular bins.
  MultiAxisSpec2D axes(std::array<AxisSpec, 2>{
      AxisSpec::Equidistant(2, 0., 2., AxisBoundaryType::Open,
                            AxisDirection::AxisX),
      AxisSpec::Equidistant(3, 0., 3., AxisBoundaryType::Open,
                            AxisDirection::AxisY)});
  auto multi = axes.buildMultiAxis();
  GridSurfaceMaterial::Direct values(20);
  for (std::size_t y = 0; y < 5; ++y) {
    for (std::size_t x = 0; x < 4; ++x) {
      values[multi->getGlobalBinFromLocalBins({x, y})] =
          MaterialSlab::Vacuum(static_cast<float>(x + 10 * y));
    }
  }
  source.surfaceMaterials.clear();
  source.surfaceMaterials[id] =
      std::make_shared<GridSurfaceMaterial>(axes, values);
  const auto j = converter.toJson(source);
  const auto& serialized = j["surfaces"][0]["material"]["storage"]["values"];
  BOOST_CHECK_EQUAL(serialized[1]["thickness"].get<double>(), 1.);
  BOOST_CHECK_EQUAL(serialized[4]["thickness"].get<double>(), 10.);
  const auto decoded = converter.fromJson(j);
  const auto& grid = dynamic_cast<const GridSurfaceMaterial&>(
      *decoded.surfaceMaterials.at(id));
  BOOST_CHECK(std::get<GridSurfaceMaterial::Direct>(grid.storage()) == values);
  BOOST_CHECK_EQUAL(grid.materialSlab(Vector2(-1., -1.)).thickness(), 0.);
  BOOST_CHECK_EQUAL(grid.materialSlab(Vector2(3., 4.)).thickness(), 43.);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentRejectsVolumes) {
  Converter converter;
  TrackingGeometryMaterial material;
  material.volumeMaterials.emplace(
      GeometryIdentifier().withVolume(1),
      std::make_shared<HomogeneousVolumeMaterial>(Material::Vacuum()));
  BOOST_CHECK_THROW(converter.toJson(material), std::invalid_argument);
  // Even an explicit null assignment must not be silently discarded.
  material.volumeMaterials.begin()->second.reset();
  BOOST_CHECK_THROW(converter.toJson(material), std::invalid_argument);
  auto document = fixture("minimal.json");
  BOOST_CHECK(!document.contains("volumes"));
  document["volumes"] = nlohmann::json::array();
  BOOST_CHECK_THROW(converter.fromJson(document), std::invalid_argument);
  document["volumes"].push_back(
      {{"geometry_id", {{"volume", 1}}}, {"material", nullptr}});
  BOOST_CHECK_THROW(converter.fromJson(document), std::invalid_argument);
  const auto decoded = converter.fromJson(fixture("minimal.json"));
  BOOST_CHECK(decoded.volumeMaterials.empty());
  BOOST_CHECK(!converter.toJson(decoded).contains("volumes"));
}

BOOST_AUTO_TEST_CASE(MaterialDocumentRejectsInvalidInputs) {
  Converter converter;
  const auto base = fixture("surfaces.json");
  auto bad = [&](auto mutation) {
    auto j = base;
    mutation(j);
    BOOST_CHECK_THROW(converter.fromJson(j), std::invalid_argument);
  };
  bad([](auto& j) { j["header"]["version"] = 2; });
  bad([](auto& j) { j["header"]["version"] = 1.5; });
  bad([](auto& j) { j["header"]["format"] = "legacy"; });
  bad([](auto& j) { j["surfaces"][0]["material"]["kind"] = "unknown"; });
  bad([](auto& j) {
    j["surfaces"][0]["target"]["geometry_id"]["volume"] = 256;
  });
  bad([](auto& j) {
    auto copy = j["surfaces"][0];
    copy["target"]["geometry_id"]["extra"] = 0;
    j["surfaces"].push_back(copy);
  });
  bad([](auto& j) {
    j["surfaces"][5]["material"]["storage"]["store"] = "missing";
  });
  bad([](auto& j) {
    j["surfaces"][4]["material"]["storage"]["indices"][4] = 2;
  });
  bad([](auto& j) {
    j["surfaces"][3]["material"]["storage"]["values"].erase(0);
  });
  bad([](auto& j) {
    j["surfaces"][3]["material"]["axes"][0]["range"] = {1., -1.};
  });
  bad([](auto& j) { j["surfaces"][5]["material"] = nullptr; });
  bad([](auto& j) {
    j["surfaces"][0]["material"]["slab"]["thickness"] =
        std::numeric_limits<double>::quiet_NaN();
  });
  auto proto = fixture("templates.json");
  proto["surfaces"][1]["material"]["material_key"] = "mismatch";
  BOOST_CHECK_THROW(converter.fromJson(proto), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentBasicParsing) {
  Converter converter;
  auto document = fixture("minimal.json");
  document.erase("header");
  BOOST_CHECK_THROW(converter.fromJson(document), std::invalid_argument);
  const nlohmann::json legacy{
      {"Surfaces", {{"entries", nlohmann::json::array()}}}};
  BOOST_CHECK_EXCEPTION(
      converter.fromJson(legacy), std::invalid_argument, [](const auto& error) {
        return std::string(error.what()).find("ActsMaterialMapMigrate") !=
               std::string::npos;
      });
  document = fixture("minimal.json");
  document["header"] = nullptr;
  BOOST_CHECK_THROW(converter.fromJson(document), nlohmann::json::exception);
  document = fixture("minimal.json");
  document["extra"] = "ignored by codec; rejected by offline schema";
  document["surfaces"][0]["material"]["extra"] = true;
  BOOST_CHECK_NO_THROW(converter.fromJson(document));
  document["surfaces"][0]["material"].erase("slab");
  BOOST_CHECK_THROW(converter.fromJson(document), nlohmann::json::exception);
  document = fixture("minimal.json");
  document["surfaces"][0]["material"]["slab"]["thickness"] = "invalid";
  BOOST_CHECK_THROW(converter.fromJson(document), nlohmann::json::exception);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentExtensionDispatchAndFiles) {
  Converter::Config config = Converter::Config::defaultConfig();
  config.encodeSurface.registerFunction(encodeCustom);
  config.decodeSurface.registerKind("application-custom", decodeCustom);
  Converter converter(std::move(config));
  TrackingGeometryMaterial source;
  source.surfaceMaterials[GeometryIdentifier().withVolume(1)] =
      std::make_shared<CustomMaterial>(7.23751);
  const auto encoded = converter.toJson(source);
  BOOST_CHECK_THROW(Converter().fromJson(encoded), std::invalid_argument);
  const auto decoded = converter.fromJson(encoded);
  BOOST_CHECK(dynamic_cast<const CustomMaterial*>(
                  decoded.surfaceMaterials.begin()->second.get()) != nullptr);
  TemporaryDirectory tmp;
  std::vector<std::string> extensions{".json", ".cbor"};
  if (detail::zstdSupported()) {
    extensions.insert(extensions.end(), {".json.zst", ".cbor.zst"});
  }
  for (const auto& ext : extensions) {
    const auto path = tmp.path() / ("material" + ext);
    converter.toFile(source, path);
    const auto renamed = tmp.path() / "unknown-extension.bin";
    std::filesystem::rename(path, renamed);
    BOOST_CHECK(converter.toJson(converter.fromFile(renamed)) == encoded);
    std::filesystem::remove(renamed);
  }
  const auto duplicate = tmp.path() / "duplicate.json";
  {
    std::ofstream out(duplicate);
    out << R"({"header":{"format":"acts-material-map","version":1,"version":1},"surfaces":[]})";
  }
  BOOST_CHECK_THROW(converter.fromFile(duplicate), std::invalid_argument);
  // CBOR map {"a":1,"a":2}: duplicate keys must also be rejected before DOM
  // construction.
  const std::vector<std::byte> cbor{
      std::byte{0xa2}, std::byte{0x61}, std::byte{0x61}, std::byte{1},
      std::byte{0x61}, std::byte{0x61}, std::byte{2}};
  BOOST_CHECK_THROW(detail::decodeJson(cbor, {}, true), std::invalid_argument);
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
