// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsFatras/Json/DetectorDescriptionJsonConverter.hpp"

#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ActsFatras::Synthetic {

namespace {

constexpr const char* kDescriptionHeader = "acts-fatras-synthetic-description";
constexpr const char* kMaterialHeader = "acts-fatras-synthetic-material";
constexpr int kFormatVersion = 0;

/// A float as the shortest decimal that reads back as itself, and an infinity
/// as the string `"inf"`, `nlohmann::json` having no number for one.
///
/// The shortest decimal matters: `nlohmann::json` holds numbers as `double`, so
/// a float widened into one prints all seventeen digits of the widening, and
/// `9.012` becomes `9.01200008392334`. These files are generated in the
/// thousands of numbers and read by people, and a table no one can scan is a
/// table no one checks.
///
/// @param value the number to write
/// @return the number, or the string "inf"
nlohmann::json encodeNumber(const float value) {
  if (std::isinf(value)) {
    return value > 0.f ? "inf" : "-inf";
  }
  if (std::isnan(value)) {
    throw std::invalid_argument(
        "DetectorDescription: a surface cannot be at a position that is not a "
        "number");
  }
  // The fewest significant digits that still read back as this float, which for
  // a float is never more than nine.
  for (int digits = 1; digits <= 9; ++digits) {
    std::array<char, 32> text{};
    std::snprintf(text.data(), text.size(), "%.*g", digits,
                  static_cast<double>(value));
    const double shortest = std::strtod(text.data(), nullptr);
    if (static_cast<float>(shortest) == value) {
      return shortest;
    }
  }
  return static_cast<double>(value);
}

/// @param j the number or the string "inf"
/// @return the number it stands for
float decodeNumber(const nlohmann::json& j) {
  if (j.is_string()) {
    const std::string text = j.get<std::string>();
    if (text == "inf") {
      return std::numeric_limits<float>::infinity();
    }
    if (text == "-inf") {
      return -std::numeric_limits<float>::infinity();
    }
    throw std::invalid_argument(
        "DetectorDescription: a bound is a number or \"inf\", not \"" + text +
        "\"");
  }
  return j.get<float>();
}

/// Look a field up, leaving the target alone where there is none.
/// @param j the object to read
/// @param key the field to look for
/// @param value where to put it
template <typename T>
void readOptional(const nlohmann::json& j, const char* key, T& value) {
  if (const auto field = j.find(key); field != j.end()) {
    value = field->get<T>();
  }
}

/// The same for a number, which may be written as the string "inf".
/// @param j the object to read
/// @param key the field to look for
/// @param value where to put it
void readOptionalNumber(const nlohmann::json& j, const char* key,
                        float& value) {
  if (const auto field = j.find(key); field != j.end()) {
    value = decodeNumber(*field);
  }
}

/// The same for a layer index, which is optional twice over: absent from the
/// file, and absent from the description because it is left to its position.
/// @param j the object to read
/// @param key the field to look for
/// @param value where to put it
void readOptional(const nlohmann::json& j, const char* key,
                  std::optional<std::uint32_t>& value) {
  if (const auto field = j.find(key); field != j.end()) {
    value = field->get<std::uint32_t>();
  }
}

/// What one band of a surface is worth in one of the six numbers a slab is made
/// of.
using BandNumber = std::function<float(const Acts::MaterialSlab&)>;

/// Write one per-band number, as a single value where the bands agree on it and
/// as one value per band otherwise.
///
/// @param j the object to write into
/// @param key the field to write
/// @param bands the bands to read
/// @param number what to read off each of them
/// @param includeVacuum whether a band of nothing has a say. It does for the two
///        lengths, a zero being how the file says a band holds nothing, and for
///        the thickness, which such a band still has. It does not for what a
///        band is made of, having nothing to be made of, and that is what keeps
///        a surface with a gap in it stated as short as it is.
void writeBandNumbers(nlohmann::json& j, const char* key,
                      const std::vector<Acts::MaterialSlab>& bands,
                      const BandNumber& number, const bool includeVacuum) {
  std::vector<float> values;
  values.reserve(bands.size());
  bool uniform = true;
  std::optional<float> shared;
  for (const Acts::MaterialSlab& band : bands) {
    const float value = number(band);
    values.push_back(value);
    if (band.isVacuum() && !includeVacuum) {
      continue;
    }
    if (!shared.has_value()) {
      shared = value;
    } else if (*shared != value) {
      uniform = false;
    }
  }
  if (uniform) {
    j[key] = encodeNumber(shared.value_or(0.f));
    return;
  }
  nlohmann::json perBand = nlohmann::json::array();
  for (const float value : values) {
    perBand.push_back(encodeNumber(value));
  }
  j[key] = std::move(perBand);
}

/// Read one per-band number, which may be stated once for the whole surface or
/// once per band.
/// @param j the object to read
/// @param key the field to read
/// @param numBands how many bands the surface has
/// @return one value per band
std::vector<float> readBandNumbers(const nlohmann::json& j, const char* key,
                                   const std::size_t numBands) {
  const nlohmann::json& field = j.at(key);
  if (!field.is_array()) {
    return std::vector<float>(numBands, decodeNumber(field));
  }
  std::vector<float> values;
  values.reserve(field.size());
  for (const nlohmann::json& value : field) {
    values.push_back(decodeNumber(value));
  }
  if (values.size() != numBands) {
    throw std::invalid_argument("SurfaceMaterial: '" + std::string(key) +
                                "' has " + std::to_string(values.size()) +
                                " values for " + std::to_string(numBands) +
                                " bands; it takes one or one per band");
  }
  return values;
}

}  // namespace

void to_json(nlohmann::json& j, const SurfaceShape& shape) {
  j = shape == SurfaceShape::Cylinder ? "cylinder" : "disc";
}

void from_json(const nlohmann::json& j, SurfaceShape& shape) {
  const std::string text = j.get<std::string>();
  if (text == "cylinder") {
    shape = SurfaceShape::Cylinder;
    return;
  }
  if (text == "disc") {
    shape = SurfaceShape::Disc;
    return;
  }
  throw std::invalid_argument(
      "DetectorDescription: a surface is a 'cylinder' or a 'disc', not '" +
      text + "'");
}

void to_json(nlohmann::json& j, const EndcapPlacement& placement) {
  j = placement == EndcapPlacement::Negative   ? "negative"
      : placement == EndcapPlacement::Positive ? "positive"
                                               : "mirrored";
}

void from_json(const nlohmann::json& j, EndcapPlacement& placement) {
  const std::string text = j.get<std::string>();
  if (text == "negative") {
    placement = EndcapPlacement::Negative;
    return;
  }
  if (text == "positive") {
    placement = EndcapPlacement::Positive;
    return;
  }
  if (text == "mirrored") {
    placement = EndcapPlacement::Mirrored;
    return;
  }
  throw std::invalid_argument(
      "DetectorDescription: an endcap is 'negative', 'positive' or 'mirrored', "
      "not '" +
      text + "'");
}

void to_json(nlohmann::json& j, const LayerKind& kind) {
  j = kind == LayerKind::Barrel   ? "barrel"
      : kind == LayerKind::Endcap ? "endcap"
                                  : "passive";
}

void from_json(const nlohmann::json& j, LayerKind& kind) {
  const std::string text = j.get<std::string>();
  if (text == "barrel") {
    kind = LayerKind::Barrel;
    return;
  }
  if (text == "endcap") {
    kind = LayerKind::Endcap;
    return;
  }
  if (text == "passive") {
    kind = LayerKind::Passive;
    return;
  }
  throw std::invalid_argument(
      "DetectorDescription: a layer is a 'barrel', an 'endcap' or a 'passive', "
      "not a '" +
      text + "'");
}

void to_json(nlohmann::json& j, const RingBounds& ring) {
  j = nlohmann::json{{"rMin", encodeNumber(ring.rMin)},
                     {"rMax", encodeNumber(ring.rMax)}};
}

void from_json(const nlohmann::json& j, RingBounds& ring) {
  ring.rMin = decodeNumber(j.at("rMin"));
  ring.rMax = decodeNumber(j.at("rMax"));
}

void to_json(nlohmann::json& j, const SurfaceMaterial& material) {
  j = nlohmann::json::object();
  if (material.bands.empty()) {
    return;
  }

  nlohmann::json bounds = nlohmann::json::array();
  for (const float edge : material.bounds) {
    bounds.push_back(encodeNumber(edge));
  }
  j["bounds"] = std::move(bounds);

  // A band of nothing is a zero in every one of these, `Acts::Material` holding
  // an infinite radiation length for one -- honest, and no number for a file.
  const auto of = [](const auto& read) {
    return [read](const Acts::MaterialSlab& band) {
      return band.isVacuum() ? 0.f : read(band.material());
    };
  };
  writeBandNumbers(j, "x0", material.bands,
                   of([](const Acts::Material& stuff) { return stuff.X0(); }),
                   true);
  writeBandNumbers(j, "l0", material.bands,
                   of([](const Acts::Material& stuff) { return stuff.L0(); }),
                   true);
  writeBandNumbers(j, "ar", material.bands,
                   of([](const Acts::Material& stuff) { return stuff.Ar(); }),
                   false);
  writeBandNumbers(j, "z", material.bands,
                   of([](const Acts::Material& stuff) { return stuff.Z(); }),
                   false);
  writeBandNumbers(j, "molarDensityX0", material.bands,
                   of([](const Acts::Material& stuff) {
                     return stuff.molarDensity() * stuff.X0();
                   }),
                   false);
  // A band of nothing still has a thickness, and it is what a crossing of it
  // costs, so it has a say in this one.
  writeBandNumbers(
      j, "thickness", material.bands,
      [](const Acts::MaterialSlab& band) { return band.thickness(); }, true);
}

void from_json(const nlohmann::json& j, SurfaceMaterial& material) {
  if (!j.contains("bounds")) {
    material = SurfaceMaterial{};
    return;
  }

  std::vector<float> bounds;
  for (const nlohmann::json& edge : j.at("bounds")) {
    bounds.push_back(decodeNumber(edge));
  }
  if (bounds.size() < 2) {
    throw std::invalid_argument(
        "SurfaceMaterial: a band needs an edge on either side of it");
  }
  const std::size_t numBands = bounds.size() - 1;

  const std::vector<float> x0 = readBandNumbers(j, "x0", numBands);
  const std::vector<float> l0 = readBandNumbers(j, "l0", numBands);
  const std::vector<float> ar = readBandNumbers(j, "ar", numBands);
  const std::vector<float> z = readBandNumbers(j, "z", numBands);
  const std::vector<float> molarDensityX0 =
      readBandNumbers(j, "molarDensityX0", numBands);
  const std::vector<float> thickness =
      readBandNumbers(j, "thickness", numBands);

  std::vector<Acts::MaterialSlab> bands;
  bands.reserve(numBands);
  for (std::size_t k = 0; k < numBands; ++k) {
    // Zero, or the infinity `Acts::Material` reports for one: either way the
    // band holds nothing and only its thickness matters.
    if (!(x0[k] > 0.f) || !(l0[k] > 0.f) || !std::isfinite(x0[k]) ||
        !std::isfinite(l0[k])) {
      bands.push_back(Acts::MaterialSlab::Vacuum(thickness[k]));
      continue;
    }
    // the density follows this band's own radiation length, which is what makes
    // one number stand for the composition of a surface banded in how much
    // material it holds
    bands.emplace_back(
        Acts::Material::fromMolarDensity(x0[k], l0[k], ar[k], z[k],
                                         molarDensityX0[k] / x0[k]),
        thickness[k]);
  }
  material = SurfaceMaterial{std::move(bounds), std::move(bands)};
}

void to_json(nlohmann::json& j, const CylinderDescription& cylinder) {
  j = nlohmann::json{
      {"radius", encodeNumber(cylinder.radius)},
      {"halfLengthZ", encodeNumber(cylinder.halfLengthZ)},
      {"modules", cylinder.modules},
      {"overlapProbability", encodeNumber(cylinder.overlapProbability)},
      {"overlapOffset", encodeNumber(cylinder.overlapOffset)}};
  if (cylinder.layer.has_value()) {
    j["layer"] = *cylinder.layer;
  }
  if (!cylinder.material.bands.empty()) {
    j["material"] = cylinder.material;
  }
}

void from_json(const nlohmann::json& j, CylinderDescription& cylinder) {
  cylinder.radius = decodeNumber(j.at("radius"));
  cylinder.halfLengthZ = decodeNumber(j.at("halfLengthZ"));
  readOptional(j, "modules", cylinder.modules);
  readOptionalNumber(j, "overlapProbability", cylinder.overlapProbability);
  readOptionalNumber(j, "overlapOffset", cylinder.overlapOffset);
  readOptional(j, "layer", cylinder.layer);
  readOptional(j, "material", cylinder.material);
}

void to_json(nlohmann::json& j, const DiscDescription& disc) {
  j = nlohmann::json{
      {"absZ", encodeNumber(disc.absZ)},
      {"rings", disc.rings},
      {"overlapProbability", encodeNumber(disc.overlapProbability)},
      {"overlapOffset", encodeNumber(disc.overlapOffset)}};
  if (disc.layer.has_value()) {
    j["layer"] = *disc.layer;
  }
  if (!disc.material.bands.empty()) {
    j["material"] = disc.material;
  }
}

void from_json(const nlohmann::json& j, DiscDescription& disc) {
  disc.absZ = decodeNumber(j.at("absZ"));
  disc.rings = j.at("rings").get<std::vector<RingBounds>>();
  readOptionalNumber(j, "overlapProbability", disc.overlapProbability);
  readOptionalNumber(j, "overlapOffset", disc.overlapOffset);
  readOptional(j, "layer", disc.layer);
  readOptional(j, "material", disc.material);
}

void to_json(nlohmann::json& j, const PassiveSurfaceDescription& passive) {
  j = nlohmann::json{{"shape", passive.shape},
                     {"refCoord", encodeNumber(passive.refCoord)},
                     {"minBound", encodeNumber(passive.minBound)},
                     {"maxBound", encodeNumber(passive.maxBound)}};
  // a cylinder straddles the interaction point and has no side to sit on
  if (passive.shape != SurfaceShape::Cylinder) {
    j["placement"] = passive.placement;
  }
  if (passive.layer.has_value()) {
    j["layer"] = *passive.layer;
  }
  if (!passive.material.bands.empty()) {
    j["material"] = passive.material;
  }
}

void from_json(const nlohmann::json& j, PassiveSurfaceDescription& passive) {
  passive.shape = j.at("shape").get<SurfaceShape>();
  passive.refCoord = decodeNumber(j.at("refCoord"));
  readOptionalNumber(j, "minBound", passive.minBound);
  readOptionalNumber(j, "maxBound", passive.maxBound);
  readOptional(j, "placement", passive.placement);
  readOptional(j, "layer", passive.layer);
  readOptional(j, "material", passive.material);
}

void to_json(nlohmann::json& j, const BarrelDescription& barrel) {
  j = nlohmann::json{{"cylinders", barrel.cylinders}};
}

void from_json(const nlohmann::json& j, BarrelDescription& barrel) {
  barrel.cylinders = j.at("cylinders").get<std::vector<CylinderDescription>>();
}

void to_json(nlohmann::json& j, const EndcapDescription& endcap) {
  j = nlohmann::json{{"placement", endcap.placement}, {"discs", endcap.discs}};
}

void from_json(const nlohmann::json& j, EndcapDescription& endcap) {
  readOptional(j, "placement", endcap.placement);
  endcap.discs = j.at("discs").get<std::vector<DiscDescription>>();
}

void to_json(nlohmann::json& j, const SubsystemDescription& subsystem) {
  j = nlohmann::json{{"name", subsystem.name}};
  if (!subsystem.barrels.empty()) {
    j["barrels"] = subsystem.barrels;
  }
  if (!subsystem.endcaps.empty()) {
    j["endcaps"] = subsystem.endcaps;
  }
  if (!subsystem.passives.empty()) {
    j["passives"] = subsystem.passives;
  }
}

void from_json(const nlohmann::json& j, SubsystemDescription& subsystem) {
  subsystem.name = j.at("name").get<std::string>();
  readOptional(j, "barrels", subsystem.barrels);
  readOptional(j, "endcaps", subsystem.endcaps);
  readOptional(j, "passives", subsystem.passives);
}

void to_json(nlohmann::json& j, const DetectorDescription& description) {
  j = nlohmann::json{{"escapeRadius", encodeNumber(description.escapeRadius)},
                     {"escapeHalfZ", encodeNumber(description.escapeHalfZ)},
                     {"subsystems", description.subsystems}};
  if (!description.passives.empty()) {
    j["passives"] = description.passives;
  }
}

void from_json(const nlohmann::json& j, DetectorDescription& description) {
  readOptionalNumber(j, "escapeRadius", description.escapeRadius);
  readOptionalNumber(j, "escapeHalfZ", description.escapeHalfZ);
  readOptional(j, "passives", description.passives);
  description.subsystems =
      j.at("subsystems").get<std::vector<SubsystemDescription>>();
}

void to_json(nlohmann::json& j, const LayerId& layer) {
  j = nlohmann::json{{"kind", layer.kind}, {"layer", layer.layer}};
  // a passive of the detector itself belongs to no subsystem
  if (!layer.subsystem.empty()) {
    j["subsystem"] = layer.subsystem;
  }
  // only an endcap is numbered per side; see `LayerId::placement`
  if (layer.kind == LayerKind::Endcap) {
    j["placement"] = layer.placement;
  }
}

void from_json(const nlohmann::json& j, LayerId& layer) {
  layer.kind = j.at("kind").get<LayerKind>();
  layer.layer = j.at("layer").get<std::uint32_t>();
  readOptional(j, "subsystem", layer.subsystem);
  readOptional(j, "placement", layer.placement);
}

void to_json(nlohmann::json& j, const MaterialEntry& entry) {
  // The identifier and the material in one object: an entry reads as "this
  // layer is made of this" rather than as a pair.
  j = entry.layer;
  j.update(nlohmann::json(entry.material));
}

void from_json(const nlohmann::json& j, MaterialEntry& entry) {
  entry.layer = j.get<LayerId>();
  entry.material = j.get<SurfaceMaterial>();
}

}  // namespace ActsFatras::Synthetic

namespace ActsFatras {

namespace {

/// @param path the file to read
/// @param header the header field it has to carry
/// @return its content
nlohmann::json readFile(const std::filesystem::path& path, const char* header) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("ActsFatras::Synthetic: cannot read " +
                             path.string());
  }
  nlohmann::json j;
  try {
    file >> j;
  } catch (const nlohmann::json::parse_error& error) {
    throw std::runtime_error("ActsFatras::Synthetic: " + path.string() +
                             " is not JSON: " + error.what());
  }
  if (!j.contains(header)) {
    throw std::runtime_error("ActsFatras::Synthetic: " + path.string() +
                             " carries no '" + header +
                             "' header, so it is some other kind of file");
  }
  return j;
}

/// @param path the file to write
/// @param header the header field to stamp it with
/// @param content what to write
void writeFile(const std::filesystem::path& path, const char* header,
               nlohmann::json content) {
  content[header] =
      nlohmann::json{{"format-version", Synthetic::kFormatVersion}};
  std::ofstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("ActsFatras::Synthetic: cannot write " +
                             path.string());
  }
  // four spaces and a single trailing newline, which is what
  // `CI/format_json.py` normalises every JSON file in the repository to
  file << content.dump(4) << std::endl;
}

}  // namespace

Synthetic::DetectorDescription Synthetic::readDetectorDescription(
    const std::filesystem::path& path) {
  return readFile(path, kDescriptionHeader).get<DetectorDescription>();
}

void Synthetic::writeDetectorDescription(
    const std::filesystem::path& path, const DetectorDescription& description) {
  writeFile(path, kDescriptionHeader, nlohmann::json(description));
}

Synthetic::MaterialDecoration Synthetic::readMaterialDecoration(
    const std::filesystem::path& path) {
  return readFile(path, kMaterialHeader).at("layers").get<MaterialDecoration>();
}

void Synthetic::writeMaterialDecoration(const std::filesystem::path& path,
                                        const MaterialDecoration& decoration) {
  writeFile(path, kMaterialHeader, nlohmann::json{{"layers", decoration}});
}

}  // namespace ActsFatras
