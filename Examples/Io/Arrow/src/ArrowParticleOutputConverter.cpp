// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowParticleOutputConverter.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

/// Nested layout: one row per event, each field a @c list<T> whose single
/// list element at row @c N holds all particles of event @c N. Uses
/// @c float32 for kinematics to match the established on-disk convention.
std::shared_ptr<arrow::Schema> particleSchema() {
  return arrow::schema({
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("pdg", arrow::list(arrow::int32()), false),
      arrow::field("charge", arrow::list(arrow::float32()), false),
      arrow::field("mass", arrow::list(arrow::float32()), false),
      arrow::field("px", arrow::list(arrow::float32()), false),
      arrow::field("py", arrow::list(arrow::float32()), false),
      arrow::field("pz", arrow::list(arrow::float32()), false),
      arrow::field("e", arrow::list(arrow::float32()), false),
      arrow::field("vx", arrow::list(arrow::float32()), false),
      arrow::field("vy", arrow::list(arrow::float32()), false),
      arrow::field("vz", arrow::list(arrow::float32()), false),
      arrow::field("vt", arrow::list(arrow::float32()), false),
  });
}

void check(const arrow::Status& s, const char* what) {
  if (!s.ok()) {
    throw std::runtime_error(std::string(what) + ": " + s.ToString());
  }
}

}  // namespace

ArrowParticleOutputConverter::ArrowParticleOutputConverter(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : ArrowOutputConverter("ArrowParticleOutputConverter", std::move(logger)),
      m_cfg(cfg) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing particles input collection");
  }
  if (m_cfg.outputTable.empty()) {
    throw std::invalid_argument("Missing output table name");
  }
  m_inputParticles.initialize(m_cfg.inputParticles);
  m_outputTable.initialize(m_cfg.outputTable);
}

std::vector<std::string> ArrowParticleOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowParticleOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const SimParticleContainer& particles = m_inputParticles(ctx);
  auto* pool = arrow::default_memory_pool();

  arrow::ListBuilder idList(pool, std::make_shared<arrow::UInt64Builder>(pool));
  arrow::ListBuilder pdgList(pool, std::make_shared<arrow::Int32Builder>(pool));
  arrow::ListBuilder chargeList(pool,
                                std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder massList(pool,
                              std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pxList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder eList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vxList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vtList(pool, std::make_shared<arrow::FloatBuilder>(pool));

  check(idList.Append(), "open particle_id list");
  check(pdgList.Append(), "open pdg list");
  check(chargeList.Append(), "open charge list");
  check(massList.Append(), "open mass list");
  check(pxList.Append(), "open px list");
  check(pyList.Append(), "open py list");
  check(pzList.Append(), "open pz list");
  check(eList.Append(), "open e list");
  check(vxList.Append(), "open vx list");
  check(vyList.Append(), "open vy list");
  check(vzList.Append(), "open vz list");
  check(vtList.Append(), "open vt list");

  auto* idV = static_cast<arrow::UInt64Builder*>(idList.value_builder());
  auto* pdgV = static_cast<arrow::Int32Builder*>(pdgList.value_builder());
  auto* chargeV =
      static_cast<arrow::FloatBuilder*>(chargeList.value_builder());
  auto* massV = static_cast<arrow::FloatBuilder*>(massList.value_builder());
  auto* pxV = static_cast<arrow::FloatBuilder*>(pxList.value_builder());
  auto* pyV = static_cast<arrow::FloatBuilder*>(pyList.value_builder());
  auto* pzV = static_cast<arrow::FloatBuilder*>(pzList.value_builder());
  auto* eV = static_cast<arrow::FloatBuilder*>(eList.value_builder());
  auto* vxV = static_cast<arrow::FloatBuilder*>(vxList.value_builder());
  auto* vyV = static_cast<arrow::FloatBuilder*>(vyList.value_builder());
  auto* vzV = static_cast<arrow::FloatBuilder*>(vzList.value_builder());
  auto* vtV = static_cast<arrow::FloatBuilder*>(vtList.value_builder());

  const auto n = particles.size();
  check(idV->Reserve(n), "reserve particle_id");
  check(pdgV->Reserve(n), "reserve pdg");
  check(chargeV->Reserve(n), "reserve charge");
  check(massV->Reserve(n), "reserve mass");
  check(pxV->Reserve(n), "reserve px");
  check(pyV->Reserve(n), "reserve py");
  check(pzV->Reserve(n), "reserve pz");
  check(eV->Reserve(n), "reserve e");
  check(vxV->Reserve(n), "reserve vx");
  check(vyV->Reserve(n), "reserve vy");
  check(vzV->Reserve(n), "reserve vz");
  check(vtV->Reserve(n), "reserve vt");

  for (const auto& particle : particles) {
    const auto& s = particle.initialState();
    const auto mom = s.momentum();
    const auto pos = s.position();
    const auto bc = s.particleId();
    const auto packedId =
        (static_cast<std::uint64_t>(bc.vertexPrimary()) << 48) |
        (static_cast<std::uint64_t>(bc.vertexSecondary()) << 32) |
        static_cast<std::uint64_t>(bc.particle());

    idV->UnsafeAppend(packedId);
    pdgV->UnsafeAppend(static_cast<std::int32_t>(s.pdg()));
    chargeV->UnsafeAppend(static_cast<float>(s.charge()));
    massV->UnsafeAppend(static_cast<float>(s.mass()));
    pxV->UnsafeAppend(static_cast<float>(mom.x()));
    pyV->UnsafeAppend(static_cast<float>(mom.y()));
    pzV->UnsafeAppend(static_cast<float>(mom.z()));
    eV->UnsafeAppend(static_cast<float>(s.energy()));
    vxV->UnsafeAppend(static_cast<float>(pos.x()));
    vyV->UnsafeAppend(static_cast<float>(pos.y()));
    vzV->UnsafeAppend(static_cast<float>(pos.z()));
    vtV->UnsafeAppend(static_cast<float>(s.time()));
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(idList),    finish(pdgList),   finish(chargeList),
      finish(massList),  finish(pxList),    finish(pyList),
      finish(pzList),    finish(eList),     finish(vxList),
      finish(vyList),    finish(vzList),    finish(vtList),
  };

  auto table = arrow::Table::Make(particleSchema(), arrays);
  m_outputTable(ctx, std::move(table));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
