// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"
#include "ActsTests/CommonHelpers/TemporaryDirectory.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
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
                            Converter::EncodeContext& context) {
  return {{"kind", "application-custom"},
          {"thickness",
           m.materialSlab(Vector2::Zero()).thickness() / context.lengthUnit()}};
}
std::unique_ptr<const ISurfaceMaterial> decodeCustom(
    const nlohmann::json& j, const Converter::DecodeContext& context) {
  return std::make_unique<CustomMaterial>(j.at("thickness").get<double>() *
                                          context.lengthUnit());
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

BOOST_AUTO_TEST_CASE(MaterialDocumentUnits) {
  Converter converter;
  auto document = fixture("minimal.json");
  document["header"]["units"] = {{"length", "cm"},
                                 {"angle", "deg"},
                                 {"energy", "MeV"},
                                 {"material_amount", "mmol"}};
  const auto originalSlab = document["surfaces"][0]["material"]["slab"];
  const auto encoded = converter.toJson(converter.fromJson(document));
  BOOST_CHECK(encoded["header"]["units"] ==
              fixture("minimal.json")["header"]["units"]);
  const auto& slab = encoded["surfaces"][0]["material"]["slab"];
  BOOST_CHECK_CLOSE_FRACTION(slab["thickness"].get<double>(),
                             originalSlab["thickness"].get<double>() * 10.,
                             1e-6);
  for (const auto* field :
       {"radiation_length", "interaction_length", "molar_density",
        "molar_electron_density", "mean_excitation_energy"}) {
    const std::string name(field);
    const double factor = name.ends_with("length")    ? 10.
                          : name.ends_with("density") ? 1e-6
                                                      : 1e-3;
    BOOST_CHECK_CLOSE_FRACTION(
        slab["material"][field].get<double>(),
        originalSlab["material"][field].get<double>() * factor, 1e-6);
  }

  // Exercise binned, direct, indexed and shared storage with noncanonical
  // units.
  auto surfaces = fixture("surfaces.json");
  const auto baseline = converter.toJson(converter.fromJson(surfaces));
  surfaces["header"]["units"] = document["header"]["units"];
  const auto rescaled = converter.toJson(converter.fromJson(surfaces));
  auto checkSlabs = [&](auto&& self, const nlohmann::json& before,
                        const nlohmann::json& after) -> void {
    if (before.is_object() && before.contains("thickness")) {
      BOOST_CHECK_CLOSE_FRACTION(after.at("thickness").get<double>(),
                                 before.at("thickness").get<double>() * 10.,
                                 1e-6);
      if (before.at("material").at("kind") == "material") {
        for (const auto* field : {"molar_density", "molar_electron_density"}) {
          BOOST_CHECK_CLOSE_FRACTION(
              after.at("material").at(field).template get<double>(),
              before.at("material").at(field).template get<double>() * 1e-6,
              1e-6);
        }
      }
    } else if (before.is_object()) {
      for (const auto& [key, value] : before.items()) {
        self(self, value, after.at(key));
      }
    } else if (before.is_array()) {
      for (std::size_t i = 0; i < before.size(); ++i) {
        self(self, before.at(i), after.at(i));
      }
    }
  };
  checkSlabs(checkSlabs, baseline, rescaled);

  auto templates = fixture("templates.json");
  templates["header"]["units"] = document["header"]["units"];
  auto& protoAxes = templates["surfaces"][1]["material"]["axes"];
  protoAxes[0] = {{"kind", "equidistant"},
                  {"bins", 2},
                  {"direction", "phi"},
                  {"range", {-180., 180.}}};
  auto& binning = templates["surfaces"][0]["material"]["binning"];
  binning["transform"] = {
      {"rotation", {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}},
      {"translation", {1., 2., 3.}}};
  const auto output = converter.toJson(converter.fromJson(templates));
  const auto& axes = output["surfaces"].back()["material"]["axes"];
  BOOST_CHECK_CLOSE_FRACTION(axes[0]["range"][1].get<double>(),
                             180. * UnitConstants::degree, 1e-12);
  BOOST_CHECK(axes[1] == protoAxes[1]);  // Normalized edges stay dimensionless.
  BOOST_CHECK(output["surfaces"][0]["material"]["binning"]["transform"]
                    ["translation"] == nlohmann::json({10., 20., 30.}));
  for (const auto* direction : {"x", "rphi", "eta", "theta"}) {
    protoAxes[0]["direction"] = direction;
    const auto converted = converter.toJson(converter.fromJson(templates));
    const double factor = std::string(direction) == "eta" ? 1.
                          : std::string(direction) == "theta"
                              ? UnitConstants::degree
                              : 10.;
    BOOST_CHECK_CLOSE_FRACTION(converted["surfaces"]
                                   .back()["material"]["axes"][0]["range"][1]
                                   .get<double>(),
                               180. * factor, 1e-12);
  }

  for (const auto* key : {"length", "angle", "energy", "material_amount"}) {
    auto invalid = document;
    invalid["header"]["units"][key] = "unknown";
    BOOST_CHECK_EXCEPTION(converter.fromJson(invalid), std::invalid_argument,
                          [&](const auto& error) {
                            return std::string(error.what()).find(key) !=
                                   std::string::npos;
                          });
    invalid["header"]["units"].erase(key);
    BOOST_CHECK_THROW(converter.fromJson(invalid), nlohmann::json::exception);
  }
  document["header"]["units"]["length"] = "GeV";
  BOOST_CHECK_THROW(converter.fromJson(document), std::invalid_argument);
  document["header"].erase("units");
  BOOST_CHECK_THROW(converter.fromJson(document), nlohmann::json::exception);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentQuantization) {
  Converter converter;
  const auto id = GeometryIdentifier().withVolume(1);
  for (unsigned int bits : {0u, 12u, 16u, 20u, 23u}) {
    Converter::Options options;
    options.materialFractionBits = bits;
    const double bound = std::ldexp(1., -static_cast<int>(bits) - 1);
    for (float value :
         {0.f, -0.f, 7.23751f, std::nextafter(1.f, 2.f),
          1.f + std::ldexp(1.f, -17), 1.f + 3 * std::ldexp(1.f, -17),
          std::numeric_limits<float>::denorm_min(),
          std::numeric_limits<float>::min(),
          std::numeric_limits<float>::max()}) {
      const auto material = Material::fromMolarDensity(
          value, value, value > 0 ? value : 28.f, value, value, value, value);
      TrackingGeometryMaterial source;
      const double split = 0.123456789012345;
      source.surfaceMaterials[id] =
          std::make_shared<HomogeneousSurfaceMaterial>(
              MaterialSlab(material, value), split);
      const auto original = converter.toJson(source);
      const auto encoded = converter.toJson(source, options);
      if (bits == 23) {
        BOOST_CHECK(encoded == original);
      }
      for (const auto& parsed :
           {nlohmann::json::parse(encoded.dump()),
            nlohmann::json::from_cbor(nlohmann::json::to_cbor(encoded))}) {
        const auto recovered = converter.fromJson(parsed);
        BOOST_CHECK(converter.toJson(recovered, options) == encoded);
        const auto& surface = *recovered.surfaceMaterials.at(id);
        const float actual = surface.materialSlab(Vector2::Zero()).thickness();
        if (!std::isnormal(value) || bits == 23) {
          BOOST_CHECK_EQUAL(std::bit_cast<std::uint32_t>(actual),
                            std::bit_cast<std::uint32_t>(value));
        } else {
          BOOST_CHECK_LE(std::abs((double(actual) - value) / value), bound);
        }
        if (bits == 16 && value == 1.f + std::ldexp(1.f, -17)) {
          BOOST_CHECK_EQUAL(actual, 1.f);
        }
        if (bits == 16 && value == 1.f + 3 * std::ldexp(1.f, -17)) {
          BOOST_CHECK_EQUAL(actual, 1.f + std::ldexp(1.f, -15));
        }
        BOOST_CHECK_EQUAL(surface.factor(Direction::Backward(),
                                         MaterialUpdateMode::PreUpdate),
                          split);
      }
      for (const auto& change : nlohmann::json::diff(original, encoded)) {
        BOOST_CHECK_EQUAL(change.at("op").get<std::string>(), "replace");
        const auto path =
            nlohmann::json::json_pointer(change.at("path").get<std::string>());
        const double before = original.at(path).get<double>();
        const double after = encoded.at(path).get<double>();
        BOOST_REQUIRE_NE(before, 0.);
        BOOST_CHECK_LE(std::abs((after - before) / before), bound);
      }
    }
  }
  Converter::Options invalid;
  invalid.materialFractionBits = 24;
  BOOST_CHECK_THROW(converter.toJson({}, invalid), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MaterialDocumentQuantizationStorageAndFiles) {
  Converter converter;
  Converter::Options options;
  options.materialFractionBits = 16;
  TemporaryDirectory tmp;
  for (const auto* name : {"minimal.json", "surfaces.json", "templates.json"}) {
    const auto source = converter.fromJson(fixture(name));
    const auto original = converter.toJson(source);
    const auto encoded = converter.toJson(source, options);
    for (const auto& change : nlohmann::json::diff(original, encoded)) {
      BOOST_CHECK_EQUAL(change.at("op").get<std::string>(), "replace");
      const auto path = change.at("path").get<std::string>();
      // IDs, indices, coordinates, settings, strings and structure must
      // survive.
      const std::array fields{"/thickness",
                              "/radiation_length",
                              "/interaction_length",
                              "/relative_atomic_mass",
                              "/atomic_number",
                              "/molar_density",
                              "/molar_electron_density",
                              "/mean_excitation_energy"};
      BOOST_CHECK(std::ranges::any_of(
          fields, [&](const auto* field) { return path.ends_with(field); }));
    }
    std::vector<std::string> extensions{".json", ".cbor"};
    if (detail::zstdSupported()) {
      extensions.insert(extensions.end(), {".json.zst", ".cbor.zst"});
    }
    for (const auto& extension : extensions) {
      const auto path = tmp.path() / (std::string(name) + extension);
      converter.toFile(source, path, options);
      BOOST_CHECK(converter.toJson(converter.fromFile(path)) == encoded);
    }
  }
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
                      .withPassive(255)
                      .withSensitive(1048575)
                      .withExtra(255);
  source.surfaceMaterials[id] = std::make_shared<HomogeneousSurfaceMaterial>(
      slab, 0.25, MappingType::Sensor);
  source.setDescription("roundtrip");
  const auto encoded = converter.toJson(source);
  BOOST_CHECK_EQUAL(encoded.at("$schema").get<std::string>(),
                    "urn:acts:material-map:1");
  const auto& encodedId =
      encoded.at("surfaces").at(0).at("target").at("geometry_id");
  BOOST_CHECK_EQUAL(encodedId.at("passive").get<unsigned int>(), 255u);
  BOOST_CHECK(!encodedId.contains("approach"));
  BOOST_CHECK_EQUAL(encodedId.at("portal").get<unsigned int>(), 255u);
  BOOST_CHECK(!encodedId.contains("boundary"));
  const auto recovered = converter.fromJson(encoded);
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
  Converter::Options quantized;
  quantized.materialFractionBits = 0;
  BOOST_CHECK(converter.toJson(source, quantized) == encoded);
  BOOST_CHECK_THROW(Converter().fromJson(encoded), std::invalid_argument);
  auto centimeters = encoded;
  centimeters["header"]["units"]["length"] = "cm";
  const auto scaled = converter.fromJson(centimeters);
  BOOST_CHECK_CLOSE_FRACTION(scaled.surfaceMaterials.begin()
                                 ->second->materialSlab(Vector2::Zero())
                                 .thickness(),
                             source.surfaceMaterials.begin()
                                     ->second->materialSlab(Vector2::Zero())
                                     .thickness() *
                                 10.,
                             1e-6);
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
}
BOOST_AUTO_TEST_SUITE_END()
}  // namespace ActsTests
