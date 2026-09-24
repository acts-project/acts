// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/TrackingGeometryMaterialJsonConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/MergedMaterialMarker.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "ActsPlugins/Json/detail/JsonIo.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>

namespace {
using namespace Acts;
using Converter = TrackingGeometryMaterialJsonConverter;
using EncodeContext = Converter::EncodeContext;
using DecodeContext = Converter::DecodeContext;
constexpr double lengthUnit = UnitConstants::mm;
constexpr double densityUnit =
    UnitConstants::mol / (lengthUnit * lengthUnit * lengthUnit);
constexpr double energyUnit = UnitConstants::GeV;

void check(bool condition, const std::string& message) {
  if (!condition) {
    throw std::invalid_argument(message);
  }
}

void array(const nlohmann::json& j, std::size_t min,
           std::size_t max = std::numeric_limits<std::size_t>::max()) {
  check(j.is_array() && j.size() >= min && j.size() <= max,
        "invalid array size");
}

float finiteFloat(double v) {
  check(std::isfinite(v) && std::abs(v) <= std::numeric_limits<float>::max(),
        "value outside finite float range");
  auto f = static_cast<float>(v);
  check(v == 0 || f != 0, "value underflows float range");
  return f;
}

std::size_t index(
    const nlohmann::json& j,
    std::size_t maximum = std::numeric_limits<std::size_t>::max()) {
  check(j.is_number_integer(), "expected integer");
  check(j.is_number_unsigned() || j.get<std::int64_t>() >= 0,
        "expected nonnegative integer");
  auto v = j.get<std::uint64_t>();
  check(v <= maximum, "integer out of range");
  return static_cast<std::size_t>(v);
}

std::size_t product(std::size_t a, std::size_t b) {
  check(b == 0 || a <= std::numeric_limits<std::size_t>::max() / b,
        "bin count overflow");
  return a * b;
}

nlohmann::json encodeId(GeometryIdentifier id) {
  nlohmann::json j = nlohmann::json::object();
  for (auto [name, value] :
       std::initializer_list<std::pair<const char*, std::uint64_t>>{
           {"volume", id.volume()},
           {"boundary", id.boundary()},
           {"layer", id.layer()},
           {"approach", id.approach()},
           {"sensitive", id.sensitive()},
           {"extra", id.extra()}}) {
    if (value != 0) {
      j[name] = value;
    }
  }
  return j;
}

GeometryIdentifier decodeId(const nlohmann::json& j) {
  auto field = [&](const char* key, std::uint64_t max) {
    return j.contains(key) ? index(j.at(key), max) : 0;
  };
  return GeometryIdentifier()
      .withVolume(field("volume", GeometryIdentifier::getMaxVolume()))
      .withBoundary(field("boundary", GeometryIdentifier::getMaxBoundary()))
      .withLayer(field("layer", GeometryIdentifier::getMaxLayer()))
      .withApproach(field("approach", GeometryIdentifier::getMaxApproach()))
      .withSensitive(field("sensitive", GeometryIdentifier::getMaxSensitive()))
      .withExtra(field("extra", GeometryIdentifier::getMaxExtra()));
}

const std::array<std::string, 9> directions{"x",    "y",     "z",   "r",  "phi",
                                            "rphi", "theta", "eta", "mag"};

AxisDirection direction(const nlohmann::json& j) {
  const auto it = std::ranges::find(directions, j.get<std::string>());
  check(it != directions.end(), "unsupported axis direction");
  return static_cast<AxisDirection>(it - directions.begin());
}

std::string direction(AxisDirection d) {
  return directions.at(static_cast<std::size_t>(d));
}

double axisUnit(std::optional<AxisDirection> d) {
  using enum AxisDirection;
  if (d == AxisPhi || d == AxisTheta || d == AxisEta) {
    return 1.;
  }
  return lengthUnit;
}

AxisBoundaryType boundary(const nlohmann::json& j) {
  using enum AxisBoundaryType;
  auto s = j.get<std::string>();
  if (s == "open") {
    return Open;
  }
  if (s == "bound") {
    return Bound;
  }
  if (s == "closed") {
    return Closed;
  }
  throw std::invalid_argument("unsupported axis boundary '" + s + "'");
}

std::string boundary(AxisBoundaryType b) {
  using enum AxisBoundaryType;
  switch (b) {
    case Open:
      return "open";
    case Bound:
      return "bound";
    case Closed:
      return "closed";
  }
  throw std::invalid_argument("invalid boundary");
}

Transform3 decodeTransform(const nlohmann::json& j) {
  array(j.at("rotation"), 3, 3);
  array(j.at("translation"), 3, 3);
  Transform3 t = Transform3::Identity();
  for (int r = 0; r < 3; ++r) {
    array(j.at("rotation").at(r), 3, 3);
    for (int c = 0; c < 3; ++c) {
      t.linear()(r, c) = j.at("rotation").at(r).at(c).get<double>();
    }
    t.translation()[r] = j.at("translation").at(r).get<double>() * lengthUnit;
  }
  return t;
}

nlohmann::json encodeTransform(const Transform3& t) {
  nlohmann::json j{{"rotation", nlohmann::json::array()},
                   {"translation", nlohmann::json::array()}};
  for (int r = 0; r < 3; ++r) {
    nlohmann::json row = nlohmann::json::array();
    for (int c = 0; c < 3; ++c) {
      row.emplace_back(t.linear()(r, c));
    }
    j["rotation"].emplace_back(row);
    j["translation"].emplace_back(t.translation()[r] / lengthUnit);
  }
  return j;
}

float decodeLength(const nlohmann::json& j) {
  if (j.is_string() && j == "infinity") {
    return std::numeric_limits<float>::infinity();
  }
  auto v = j.get<double>();
  return finiteFloat(v * lengthUnit);
}

nlohmann::json encodeLength(float v) {
  if (v == std::numeric_limits<float>::infinity()) {
    return "infinity";
  }
  return v / lengthUnit;
}

Material decodeMaterial(const nlohmann::json& j) {
  auto kind = j.at("kind").get<std::string>();
  if (kind == "vacuum") {
    return Material::Vacuum();
  }
  check(kind == "material", "unsupported composition kind");
  auto ar = finiteFloat(j.at("relative_atomic_mass").get<double>());
  return Material::fromMolarDensity(
      decodeLength(j.at("radiation_length")),
      decodeLength(j.at("interaction_length")), ar,
      finiteFloat(j.at("atomic_number").get<double>()),
      finiteFloat(j.at("molar_density").get<double>() * densityUnit),
      finiteFloat(j.at("molar_electron_density").get<double>() * densityUnit),
      finiteFloat(j.at("mean_excitation_energy").get<double>() * energyUnit));
}

nlohmann::json encodeMaterial(const Material& m) {
  if (m.isVacuum()) {
    return nlohmann::json{{"kind", "vacuum"}};
  }
  nlohmann::json j{
      {"kind", "material"},
      {"radiation_length", encodeLength(m.X0())},
      {"interaction_length", encodeLength(m.L0())},
      {"relative_atomic_mass", m.Ar()},
      {"atomic_number", m.Z()},
      {"molar_density", m.molarDensity() / densityUnit},
      {"molar_electron_density", m.molarElectronDensity() / densityUnit},
      {"mean_excitation_energy", m.meanExcitationEnergy() / energyUnit}};
  return j;
}

MaterialSlab decodeSlab(const nlohmann::json& j) {
  return MaterialSlab(
      decodeMaterial(j.at("material")),
      finiteFloat(j.at("thickness").get<double>() * lengthUnit));
}

nlohmann::json encodeSlab(const MaterialSlab& slab) {
  return {{"material", encodeMaterial(slab.material())},
          {"thickness", slab.thickness() / lengthUnit}};
}

std::vector<MaterialSlab> decodeSlabs(const nlohmann::json& j) {
  array(j, 1);
  std::vector<MaterialSlab> result;
  result.reserve(j.size());
  for (const auto& v : j) {
    result.emplace_back(decodeSlab(v));
  }
  return result;
}

nlohmann::json encodeSlabs(const std::vector<MaterialSlab>& slabs) {
  check(!slabs.empty(), "empty slab store");
  nlohmann::json result = nlohmann::json::array();
  for (const auto& s : slabs) {
    result.emplace_back(encodeSlab(s));
  }
  return result;
}

constexpr std::array<std::string_view, 4> mappingNames{"pre", "default", "post",
                                                       "sensor"};

std::pair<double, MappingType> settings(const nlohmann::json& j) {
  const auto name = j.at("mapping_type").get<std::string>();
  const auto it = std::ranges::find(mappingNames, name);
  check(it != mappingNames.end(), "unsupported mapping type");
  return {j.at("split_factor").get<double>(),
          static_cast<MappingType>(it - mappingNames.begin() - 1)};
}

nlohmann::json settings(const ISurfaceMaterial& m) {
  return {
      {"mapping_type", mappingNames.at(static_cast<int>(m.mappingType()) + 1)},
      {"split_factor",
       m.factor(Direction::Backward(), MaterialUpdateMode::PreUpdate)}};
}

AxisSpec decodeAxis(const nlohmann::json& j, bool deferred) {
  const auto kind = j.at("kind").get<std::string>();
  std::optional<AxisDirection> d;
  std::optional<AxisBoundaryType> b;
  if (j.contains("direction")) {
    d = direction(j.at("direction"));
  }
  if (j.contains("boundary")) {
    b = boundary(j.at("boundary"));
  }
  check(deferred || b.has_value(), "resolved axis requires boundary");
  if (kind == "equidistant") {
    const auto n =
        index(j.at("bins"), std::numeric_limits<std::size_t>::max() - 2);
    std::optional<double> min;
    std::optional<double> max;
    if (j.contains("range")) {
      array(j.at("range"), 2, 2);
      min = j.at("range").at(0).get<double>() * axisUnit(d);
      max = j.at("range").at(1).get<double>() * axisUnit(d);
    }
    check(deferred || min.has_value(), "resolved axis requires range");
    return AxisSpec::Equidistant(n, min, max, b, d);
  }
  check(kind == "variable" || (deferred && kind == "deferred-variable"),
        "unsupported axis kind");
  const std::string field = kind == "variable" ? "edges" : "normalized_edges";
  array(j.at(field), 2);
  std::vector<double> edges;
  for (const auto& e : j.at(field)) {
    edges.emplace_back(e.get<double>() *
                       (kind == "variable" ? axisUnit(d) : 1.));
  }
  return kind == "variable" ? AxisSpec::Variable(edges, b, d)
                            : AxisSpec::DeferredVariable(edges, b, d);
}

nlohmann::json encodeAxis(const AxisSpec& axis) {
  nlohmann::json j;
  if (axis.isEquidistant()) {
    const auto& a = axis.asEquidistant();
    check(a.min.has_value() == a.max.has_value(),
          "partially specified range is not representable");
    j = {{"kind", "equidistant"}, {"bins", a.nBins}};
    if (a.min.has_value()) {
      j["range"] = {*a.min / axisUnit(axis.direction()),
                    *a.max / axisUnit(axis.direction())};
    }
  } else if (axis.isDeferredVariable()) {
    j = {{"kind", "deferred-variable"},
         {"normalized_edges", axis.asDeferredVariable().normalizedEdges}};
  } else {
    auto edges = axis.asVariable().edges;
    for (auto& v : edges) {
      v /= axisUnit(axis.direction());
    }
    j = {{"kind", "variable"}, {"edges", edges}};
  }
  if (const auto d = axis.direction(); d) {
    j["direction"] = direction(*d);
  }
  if (const auto b = axis.boundaryType(); b) {
    j["boundary"] = boundary(*b);
  }
  return j;
}

BinningData decodeBinAxis(const nlohmann::json& j, unsigned int depth = 0) {
  using enum AxisDirection;
  check(depth < 32, "binning refinement nesting exceeds 32");
  if (j.at("kind") != "subdivided") {
    auto spec = decodeAxis(j, false);
    check(spec.direction().has_value(), "BinUtility axis requires direction");
    check(spec.boundaryType() != AxisBoundaryType::Open,
          "BinUtility has no guard cells");
    // Theta and magnitude have historical global-projection semantics that do
    // not match the physical coordinate names in the new format.
    check(spec.direction() != AxisTheta && spec.direction() != AxisMag,
          "legacy theta/mag BinUtility projection is not supported by this "
          "format");
    auto axis = spec.buildAxis();
    std::optional<float> previous;
    for (double e : axis->getBinEdges()) {
      const float edge = finiteFloat(e);
      check(!previous || edge > *previous,
            "BinUtility edges collapse at float precision");
      previous = edge;
    }
    return BinningData(*axis);
  }
  check(j.at("base").at("kind") != "subdivided",
        "refinement base must be resolved");
  auto base = decodeBinAxis(j.at("base"), depth + 1);
  auto sub = decodeBinAxis(j.at("subdivision"), depth + 1);
  check(base.binvalue == sub.binvalue && base.option == sub.option,
        "refinement axis mismatch");
  auto mode = j.at("mode").get<std::string>();
  check(mode == "replace" || mode == "repeat", "invalid refinement mode");
  const auto& edges = base.boundaries();
  if (mode == "replace") {
    bool match = false;
    for (std::size_t i = 1; i < edges.size(); ++i) {
      match |= sub.min == edges[i - 1] && sub.max == edges[i];
    }
    check(match, "replacement must match one base interval");
    check(base.bins() <= std::numeric_limits<std::size_t>::max() - sub.bins(),
          "refinement count overflow");
  } else {
    check(base.type == equidistant && sub.min == base.min &&
              sub.max == edges.at(1),
          "repeat needs an equidistant base and first-interval subdivision");
    product(base.bins(), sub.bins());
  }
  auto child = std::make_unique<const BinningData>(sub);
  if (base.type == equidistant) {
    return BinningData(base.option, base.binvalue, base.bins(), base.min,
                       base.max, std::move(child), mode == "replace");
  }
  return BinningData(base.option, base.binvalue, edges, std::move(child));
}

nlohmann::json encodeBinAxis(const BinningData& b) {
  BinningData base(b);
  base.subBinningData.reset();
  const auto& rawEdges = base.boundaries();
  nlohmann::json j{{"boundary", b.option == closed ? "closed" : "bound"},
                   {"direction", direction(b.binvalue)}};
  if (b.type == equidistant) {
    j["kind"] = "equidistant";
    j["bins"] = rawEdges.size() - 1;
    j["range"] = {b.min / axisUnit(b.binvalue), b.max / axisUnit(b.binvalue)};
  } else {
    j["kind"] = "variable";
    j["edges"] = nlohmann::json::array();
    for (float e : rawEdges) {
      j["edges"].emplace_back(e / axisUnit(b.binvalue));
    }
  }
  if (b.subBinningData) {
    j = {{"kind", "subdivided"},
         {"base", j},
         {"mode", b.subBinningAdditive ? "replace" : "repeat"},
         {"subdivision", encodeBinAxis(*b.subBinningData)}};
  }
  const auto decoded = decodeBinAxis(j);
  check(decoded.bins() == b.bins(), "refinement bin count cannot be preserved");
  return j;
}

BinUtility decodeBinning(const nlohmann::json& j, std::size_t min,
                         std::size_t max) {
  array(j.at("axes"), min, max);
  BinUtility b(j.contains("transform") ? decodeTransform(j.at("transform"))
                                       : Transform3::Identity());
  for (const auto& a : j.at("axes")) {
    b += BinUtility(decodeBinAxis(a));
  }
  return b;
}

nlohmann::json encodeBinning(const BinUtility& b) {
  nlohmann::json j{{"axes", nlohmann::json::array()}};
  for (const auto& a : b.binningData()) {
    j["axes"].emplace_back(encodeBinAxis(a));
  }
  if (!b.transform().isApprox(Transform3::Identity())) {
    j["transform"] = encodeTransform(b.transform());
  }
  return j;
}

std::optional<std::string> materialKey(const nlohmann::json& j) {
  return j.contains("material_key")
             ? std::optional(j.at("material_key").get<std::string>())
             : std::nullopt;
}

void checkKey(const ISurfaceMaterial& material, const std::string& key) {
  auto checkProto = [&](const auto* p) {
    if (p != nullptr && p->materialKey()) {
      check(*p->materialKey() == key,
            "proto material key disagrees with assignment");
    }
  };
  checkProto(dynamic_cast<const ProtoSurfaceMaterial*>(&material));
  checkProto(dynamic_cast<const ProtoGridSurfaceMaterial*>(&material));
}

nlohmann::json encodeHomogeneousSurface(const HomogeneousSurfaceMaterial& m,
                                        EncodeContext& /*context*/) {
  return {{"kind", "homogeneous"},
          {"settings", settings(m)},
          {"slab", encodeSlab(m.materialSlab())}};
}

std::unique_ptr<const ISurfaceMaterial> decodeHomogeneousSurface(
    const nlohmann::json& j, const DecodeContext& /*context*/) {
  auto [split, mapping] = settings(j.at("settings"));
  return std::make_unique<HomogeneousSurfaceMaterial>(decodeSlab(j.at("slab")),
                                                      split, mapping);
}

nlohmann::json encodeBinned(const BinnedSurfaceMaterial& m,
                            EncodeContext& /*context*/) {
  nlohmann::json values = nlohmann::json::array();
  for (const auto& row : m.fullMaterial()) {
    for (const auto& s : row) {
      values.emplace_back(encodeSlab(s));
    }
  }
  return {{"kind", "binned"},
          {"settings", settings(m)},
          {"binning", encodeBinning(m.binUtility())},
          {"values", values}};
}

std::unique_ptr<const ISurfaceMaterial> decodeBinned(
    const nlohmann::json& j, const DecodeContext& /*context*/) {
  auto [split, mapping] = settings(j.at("settings"));
  auto bins = decodeBinning(j.at("binning"), 1, 2);
  auto n0 = bins.binningData()[0].bins();
  auto n1 = bins.dimensions() == 2 ? bins.binningData()[1].bins() : 1;
  auto count = product(n0, n1);
  array(j.at("values"), count, count);
  MaterialSlabMatrix matrix(n1);
  for (std::size_t i1 = 0; i1 < n1; ++i1) {
    for (std::size_t i0 = 0; i0 < n0; ++i0) {
      matrix[i1].emplace_back(decodeSlab(j.at("values").at(i0 + n0 * i1)));
    }
  }
  return std::make_unique<BinnedSurfaceMaterial>(bins, std::move(matrix), split,
                                                 mapping);
}

nlohmann::json encodeProtoSurface(const ProtoSurfaceMaterial& m,
                                  EncodeContext& /*context*/) {
  nlohmann::json j{{"kind", "proto"},
                   {"settings", settings(m)},
                   {"binning", encodeBinning(m.binning())}};
  if (const auto& key = m.materialKey(); key) {
    j["material_key"] = *key;
  }
  return j;
}

std::unique_ptr<const ISurfaceMaterial> decodeProtoSurface(
    const nlohmann::json& j, const DecodeContext& /*context*/) {
  auto [split, mapping] = settings(j.at("settings"));
  check(split == 1, "proto surface split factor must be one");
  return std::make_unique<ProtoSurfaceMaterial>(
      decodeBinning(j.at("binning"), 0, 2), mapping, materialKey(j));
}

nlohmann::json encodeProtoGrid(const ProtoGridSurfaceMaterial& m,
                               EncodeContext& /*context*/) {
  nlohmann::json j{{"kind", "proto-grid"},
                   {"settings", settings(m)},
                   {"axes", nlohmann::json::array()}};
  for (const auto& a : m.binning().axisSpecs()) {
    j["axes"].emplace_back(encodeAxis(a));
  }
  if (const auto& key = m.materialKey(); key) {
    j["material_key"] = *key;
  }
  return j;
}

std::unique_ptr<const ISurfaceMaterial> decodeProtoGrid(
    const nlohmann::json& j, const DecodeContext& /*context*/) {
  array(j.at("axes"), 2, 2);
  auto [split, mapping] = settings(j.at("settings"));
  check(split == 1, "proto grid split factor must be one");
  MultiAxisSpec2D axes(
      std::array<AxisSpec, 2>{decodeAxis(j.at("axes").at(0), true),
                              decodeAxis(j.at("axes").at(1), true)});
  return std::make_unique<ProtoGridSurfaceMaterial>(axes, mapping,
                                                    materialKey(j));
}

nlohmann::json encodeMarker(const MergedMaterialMarker& m,
                            EncodeContext& /*context*/) {
  nlohmann::json origins = nlohmann::json::array();
  for (const auto& origin : m.origins()) {
    nlohmann::json j{{"geometry_id", encodeId(origin.geometryId)}};
    if (origin.materialKey) {
      check(!origin.materialKey->empty(), "empty origin key");
      j["material_key"] = *origin.materialKey;
    }
    origins.emplace_back(j);
  }
  return {{"kind", "merged-material-marker"}, {"origins", origins}};
}

std::unique_ptr<const ISurfaceMaterial> decodeMarker(
    const nlohmann::json& j, const DecodeContext& /*context*/) {
  array(j.at("origins"), 0);
  std::vector<MergedMaterialMarker::Origin> origins;
  for (const auto& origin : j.at("origins")) {
    origins.emplace_back(decodeId(origin.at("geometry_id")),
                         materialKey(origin));
  }
  return std::make_unique<MergedMaterialMarker>(origins);
}

nlohmann::json encodeGrid(const GridSurfaceMaterial& m,
                          EncodeContext& context) {
  nlohmann::json axes = nlohmann::json::array();
  for (const auto& a : m.binning().axisSpecs()) {
    axes.emplace_back(encodeAxis(a));
  }
  nlohmann::json storage;
  const auto n = m.multiAxis().getNBins();
  std::visit(
      [&]<typename T>(const T& s) {
        using Storage = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<Storage, GridSurfaceMaterial::Direct>) {
          storage = {{"kind", "direct"}, {"values", nlohmann::json::array()}};
        } else if constexpr (std::is_same_v<Storage,
                                            GridSurfaceMaterial::Indexed>) {
          storage = {{"kind", "indexed"},
                     {"slabs", encodeSlabs(s.material)},
                     {"indices", nlohmann::json::array()}};
        } else {
          storage = {{"kind", "globally-indexed"},
                     {"store", context.storeId(s.material)},
                     {"indices", nlohmann::json::array()}};
        }
        for (std::size_t i1 = 0; i1 < n[1] + 2; ++i1) {
          for (std::size_t i0 = 0; i0 < n[0] + 2; ++i0) {
            auto bin = m.multiAxis().getGlobalBinFromLocalBins({i0, i1});
            if constexpr (std::is_same_v<Storage,
                                         GridSurfaceMaterial::Direct>) {
              storage["values"].emplace_back(encodeSlab(s.at(bin)));
            } else {
              const auto size = [&] {
                if constexpr (std::is_same_v<Storage,
                                             GridSurfaceMaterial::Indexed>) {
                  return s.material.size();
                } else {
                  return s.material->size();
                }
              }();
              check(s.indices.at(bin) < size, "grid slab index out of range");
              storage["indices"].emplace_back(s.indices.at(bin));
            }
          }
        }
      },
      m.storage());
  return {{"kind", "grid"},
          {"settings", settings(m)},
          {"axes", axes},
          {"storage", storage}};
}

std::unique_ptr<const ISurfaceMaterial> decodeGrid(
    const nlohmann::json& j, const DecodeContext& context) {
  array(j.at("axes"), 2, 2);
  auto [split, mapping] = settings(j.at("settings"));
  MultiAxisSpec2D spec(
      std::array<AxisSpec, 2>{decodeAxis(j.at("axes").at(0), false),
                              decodeAxis(j.at("axes").at(1), false)});
  auto axes = spec.buildMultiAxis();
  const auto n = axes->getNBins();
  auto count = product(n[0] + 2, n[1] + 2);
  const auto& storage = j.at("storage");
  GridSurfaceMaterial::Storage result;
  if (const auto kind = storage.at("kind").get<std::string>();
      kind == "direct") {
    array(storage.at("values"), count, count);
    GridSurfaceMaterial::Direct slabs(count);
    for (std::size_t i1 = 0; i1 < n[1] + 2; ++i1) {
      for (std::size_t i0 = 0; i0 < n[0] + 2; ++i0) {
        slabs[axes->getGlobalBinFromLocalBins({i0, i1})] =
            decodeSlab(storage.at("values").at(i0 + (n[0] + 2) * i1));
      }
    }
    result = std::move(slabs);
  } else {
    Converter::SlabStore shared;
    std::vector<MaterialSlab> local;
    if (kind == "indexed") {
      local = decodeSlabs(storage.at("slabs"));
    } else {
      check(kind == "globally-indexed", "unknown grid storage kind");
      const auto name = storage.at("store").get<std::string>();
      shared = context.store(name);
    }
    const auto size = shared ? shared->size() : local.size();
    check(size != 0, "empty slab store");
    array(storage.at("indices"), count, count);
    std::vector<std::size_t> indices(count);
    for (std::size_t i1 = 0; i1 < n[1] + 2; ++i1) {
      for (std::size_t i0 = 0; i0 < n[0] + 2; ++i0) {
        indices[axes->getGlobalBinFromLocalBins({i0, i1})] =
            index(storage.at("indices").at(i0 + (n[0] + 2) * i1), size - 1);
      }
    }
    if (shared) {
      result = GridSurfaceMaterial::GloballyIndexed{std::move(indices),
                                                    std::move(shared)};
    } else {
      result =
          GridSurfaceMaterial::Indexed{std::move(indices), std::move(local)};
    }
  }
  return std::make_unique<GridSurfaceMaterial>(
      std::move(spec), std::move(result), split, mapping);
}

}  // namespace

namespace Acts {

std::string TrackingGeometryMaterialJsonConverter::EncodeContext::storeId(
    const SlabStore& store) {
  check(store != nullptr && !store->empty(), "null or empty shared slab store");
  for (const auto& [id, existing] : m_stores) {
    if (existing.get() == store.get()) {
      return id;
    }
  }
  std::string id = "store-" + std::to_string(m_stores.size());
  while (m_stores.contains(id)) {
    id += "-";
  }
  m_stores.try_emplace(id, store);
  return id;
}

TrackingGeometryMaterialJsonConverter::SlabStore
TrackingGeometryMaterialJsonConverter::DecodeContext::store(
    const std::string& name) const {
  const auto found = m_stores.find(name);
  check(found != m_stores.end(), "unresolved slab store '" + name + "'");
  return found->second;
}

TrackingGeometryMaterialJsonConverter::Config
TrackingGeometryMaterialJsonConverter::Config::defaultConfig() {
  Config c;
  c.encodeSurface.registerFunction(encodeHomogeneousSurface)
      .registerFunction(encodeBinned)
      .registerFunction(encodeGrid)
      .registerFunction(encodeProtoSurface)
      .registerFunction(encodeProtoGrid)
      .registerFunction(encodeMarker);
  c.decodeSurface.registerKind("homogeneous", decodeHomogeneousSurface)
      .registerKind("binned", decodeBinned)
      .registerKind("grid", decodeGrid)
      .registerKind("proto", decodeProtoSurface)
      .registerKind("proto-grid", decodeProtoGrid)
      .registerKind("merged-material-marker", decodeMarker);
  return c;
}

TrackingGeometryMaterialJsonConverter::TrackingGeometryMaterialJsonConverter(
    Config config)
    : m_config(std::move(config)) {}

nlohmann::json TrackingGeometryMaterialJsonConverter::toJson(
    const TrackingGeometryMaterial& material) const {
  check(material.volumeMaterials.empty(),
        "material document version 1 supports surface material only; "
        "volume assignments cannot be serialized");
  nlohmann::json j{
      {"$schema", "urn:acts:material-map:1-draft"},
      {"header", {{"format", "acts-material-map"}, {"version", 1}}},
      {"surfaces", nlohmann::json::array()}};
  if (material.description()) {
    j["header"]["description"] = *material.description();
  }
  EncodeContext context;
  for (const auto& [id, payload] : material.surfaceMaterials) {
    nlohmann::json value = payload ? m_config.encodeSurface(*payload, context)
                                   : nlohmann::json(nullptr);
    j["surfaces"].emplace_back(nlohmann::json{
        {"target", {{"kind", "geometry-id"}, {"geometry_id", encodeId(id)}}},
        {"material", value}});
  }
  for (const auto& [key, assignment] : material.keyedSurfaces) {
    check(!key.empty() && assignment.material != nullptr,
          "keyed assignment needs nonempty key and material");
    checkKey(*assignment.material, key);
    auto value = m_config.encodeSurface(*assignment.material, context);
    j["surfaces"].emplace_back(nlohmann::json{
        {"target",
         {{"kind", "stable-key"},
          {"key", key},
          {"recorded_geometry_id", encodeId(assignment.geometryId)}}},
        {"material", value}});
  }
  if (!context.m_stores.empty()) {
    j["slab_stores"] = nlohmann::json::object();
    for (const auto& [id, store] : context.m_stores) {
      check(store != nullptr, "null slab store");
      j["slab_stores"][id] = encodeSlabs(*store);
    }
  }
  return j;
}

TrackingGeometryMaterial TrackingGeometryMaterialJsonConverter::fromJson(
    const nlohmann::json& encoded) const {
  check(!encoded.contains("volumes"), "version 1 does not support volumes");
  if (encoded.contains("Surfaces") || encoded.contains("Volumes") ||
      encoded.contains("acts-geometry-hierarchy-map")) {
    throw std::invalid_argument(
        "Legacy material format is not supported by "
        "TrackingGeometryMaterialJsonConverter. "
        "Convert the file with ActsMaterialMapMigrate <input> <output> first.");
  }
  check(
      encoded.contains("header"),
      "Material document is missing the required header (format and version)");
  const auto& header = encoded.at("header");
  check(header.at("format") == "acts-material-map",
        "unsupported material format");
  check(index(header.at("version")) == 1,
        "unsupported material document version");
  TrackingGeometryMaterial result;
  if (header.contains("description")) {
    result.setDescription(header.at("description").get<std::string>());
  }
  DecodeContext context;
  if (encoded.contains("slab_stores")) {
    for (const auto& [name, slabs] :
         encoded.at("slab_stores").get_ref<const nlohmann::json::object_t&>()) {
      check(!name.empty(), "empty slab store name");
      context.m_stores.try_emplace(
          name,
          std::make_shared<std::vector<MaterialSlab>>(decodeSlabs(slabs)));
    }
  }
  array(encoded.at("surfaces"), 0);
  for (const auto& entry : encoded.at("surfaces")) {
    const auto& target = entry.at("target");
    const auto kind = target.at("kind").get<std::string>();
    std::shared_ptr<const ISurfaceMaterial> payload;
    if (!entry.at("material").is_null()) {
      payload = m_config.decodeSurface(entry.at("material"), context);
      check(payload != nullptr,
            "decoder returned null for a non-null surface payload");
    }
    if (kind == "geometry-id") {
      check(result.surfaceMaterials
                .try_emplace(decodeId(target.at("geometry_id")),
                             std::move(payload))
                .second,
            "duplicate surface geometry ID");
    } else {
      check(kind == "stable-key", "unsupported surface target kind");
      const auto key = target.at("key").get<std::string>();
      check(payload != nullptr, "keyed material must not be null");
      checkKey(*payload, key);
      check(result.keyedSurfaces
                .try_emplace(key, decodeId(target.at("recorded_geometry_id")),
                             std::move(payload))
                .second,
            "duplicate stable key");
    }
  }
  return result;
}

void TrackingGeometryMaterialJsonConverter::toFile(
    const TrackingGeometryMaterial& material, const std::filesystem::path& path,
    const Options& options) const {
  detail::writeJsonFile(path, toJson(material), options.indentation,
                        options.compressionLevel);
}

TrackingGeometryMaterial TrackingGeometryMaterialJsonConverter::fromFile(
    const std::filesystem::path& path) const {
  return fromJson(detail::readJsonFile(path, true));
}
}  // namespace Acts
