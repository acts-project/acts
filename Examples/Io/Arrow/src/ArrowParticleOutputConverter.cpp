// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowParticleOutputConverter.hpp"

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

namespace {

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

ArrowParticleOutputConverter::~ArrowParticleOutputConverter() = default;

std::vector<std::string> ArrowParticleOutputConverter::collections() const {
  return {m_cfg.outputTable};
}

ProcessCode ArrowParticleOutputConverter::execute(
    const AlgorithmContext& ctx) const {
  const SimParticleContainer& particles = m_inputParticles(ctx);
  auto* pool = arrow::default_memory_pool();

  arrow::ListBuilder idList(pool, std::make_shared<arrow::UInt64Builder>(pool));
  arrow::ListBuilder pdgList(pool, std::make_shared<arrow::Int64Builder>(pool));
  arrow::ListBuilder massList(pool,
                              std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder energyList(pool,
                                std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder chargeList(pool,
                                std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vxList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder timeList(pool,
                              std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pxList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pyList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder pzList(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder d0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder z0List(pool, std::make_shared<arrow::FloatBuilder>(pool));
  arrow::ListBuilder vprimList(pool,
                               std::make_shared<arrow::UInt16Builder>(pool));
  arrow::ListBuilder parentList(pool,
                                std::make_shared<arrow::Int64Builder>(pool));
  arrow::ListBuilder primaryList(pool,
                                 std::make_shared<arrow::BooleanBuilder>(pool));

  check(idList.Append(), "open particle_id list");
  check(pdgList.Append(), "open pdg_id list");
  check(massList.Append(), "open mass list");
  check(energyList.Append(), "open energy list");
  check(chargeList.Append(), "open charge list");
  check(vxList.Append(), "open vx list");
  check(vyList.Append(), "open vy list");
  check(vzList.Append(), "open vz list");
  check(timeList.Append(), "open time list");
  check(pxList.Append(), "open px list");
  check(pyList.Append(), "open py list");
  check(pzList.Append(), "open pz list");
  check(d0List.Append(), "open perigee_d0 list");
  check(z0List.Append(), "open perigee_z0 list");
  check(vprimList.Append(), "open vertex_primary list");
  check(parentList.Append(), "open parent_id list");
  check(primaryList.Append(), "open primary list");

  auto* idV = static_cast<arrow::UInt64Builder*>(idList.value_builder());
  auto* pdgV = static_cast<arrow::Int64Builder*>(pdgList.value_builder());
  auto* massV = static_cast<arrow::FloatBuilder*>(massList.value_builder());
  auto* energyV = static_cast<arrow::FloatBuilder*>(energyList.value_builder());
  auto* chargeV = static_cast<arrow::FloatBuilder*>(chargeList.value_builder());
  auto* vxV = static_cast<arrow::FloatBuilder*>(vxList.value_builder());
  auto* vyV = static_cast<arrow::FloatBuilder*>(vyList.value_builder());
  auto* vzV = static_cast<arrow::FloatBuilder*>(vzList.value_builder());
  auto* timeV = static_cast<arrow::FloatBuilder*>(timeList.value_builder());
  auto* pxV = static_cast<arrow::FloatBuilder*>(pxList.value_builder());
  auto* pyV = static_cast<arrow::FloatBuilder*>(pyList.value_builder());
  auto* pzV = static_cast<arrow::FloatBuilder*>(pzList.value_builder());
  auto* d0V = static_cast<arrow::FloatBuilder*>(d0List.value_builder());
  auto* z0V = static_cast<arrow::FloatBuilder*>(z0List.value_builder());
  auto* vprimV = static_cast<arrow::UInt16Builder*>(vprimList.value_builder());
  auto* parentV = static_cast<arrow::Int64Builder*>(parentList.value_builder());
  auto* primaryV =
      static_cast<arrow::BooleanBuilder*>(primaryList.value_builder());

  const auto n = particles.size();
  check(idV->Reserve(n), "reserve particle_id");
  check(pdgV->Reserve(n), "reserve pdg_id");
  check(massV->Reserve(n), "reserve mass");
  check(energyV->Reserve(n), "reserve energy");
  check(chargeV->Reserve(n), "reserve charge");
  check(vxV->Reserve(n), "reserve vx");
  check(vyV->Reserve(n), "reserve vy");
  check(vzV->Reserve(n), "reserve vz");
  check(timeV->Reserve(n), "reserve time");
  check(pxV->Reserve(n), "reserve px");
  check(pyV->Reserve(n), "reserve py");
  check(pzV->Reserve(n), "reserve pz");
  check(d0V->Reserve(n), "reserve perigee_d0");
  check(z0V->Reserve(n), "reserve perigee_z0");
  check(vprimV->Reserve(n), "reserve vertex_primary");
  check(parentV->Reserve(n), "reserve parent_id");
  check(primaryV->Reserve(n), "reserve primary");

  std::int64_t rowIndex = 0;
  for (const auto& particle : particles) {
    const auto& s = particle.initialState();
    const auto mom = s.momentum();
    const auto pos = s.position();
    const auto bc = s.particleId();

    // Emit the row index as the particle id (matches the colliderml
    // convention). Indices are stable within this event/output table even
    // when upstream filtering has dropped some EDM4hep particles.
    idV->UnsafeAppend(static_cast<std::uint64_t>(rowIndex));
    pdgV->UnsafeAppend(static_cast<std::int64_t>(s.pdg()));
    massV->UnsafeAppend(static_cast<float>(s.mass()));
    energyV->UnsafeAppend(static_cast<float>(s.energy()));
    chargeV->UnsafeAppend(static_cast<float>(s.charge()));
    vxV->UnsafeAppend(static_cast<float>(pos.x()));
    vyV->UnsafeAppend(static_cast<float>(pos.y()));
    vzV->UnsafeAppend(static_cast<float>(pos.z()));
    timeV->UnsafeAppend(static_cast<float>(s.time()));
    pxV->UnsafeAppend(static_cast<float>(mom.x()));
    pyV->UnsafeAppend(static_cast<float>(mom.y()));
    pzV->UnsafeAppend(static_cast<float>(mom.z()));

    // Perigee parameters are not computed here yet: the truth-to-perigee
    // propagation will be added back in a follow-up PR. Until then these
    // columns are emitted as null for every particle.
    d0V->UnsafeAppendNull();
    z0V->UnsafeAppendNull();

    vprimV->UnsafeAppend(static_cast<std::uint16_t>(bc.vertexPrimary()));
    // Emit the parent's row index in this same table so consumers can walk
    // the chain. -1 means "unknown" (parent was filtered out, or simulation
    // engine didn't record it). The container is a flat_set sorted by
    // barcode, so find()+distance is O(log N) with O(1) random access.
    const auto parentBc = particle.parentParticleId();
    std::int64_t parentRow = -1;
    if (parentBc.isValid()) {
      auto it = particles.find(parentBc);
      if (it != particles.end()) {
        parentRow = std::distance(particles.begin(), it);
      }
    }
    parentV->UnsafeAppend(parentRow);
    check(primaryV->Append(bc.generation() == 0), "append primary");
    ++rowIndex;
  }

  auto finish = [](arrow::ListBuilder& b) {
    std::shared_ptr<arrow::Array> out;
    check(b.Finish(&out), "finish list");
    return out;
  };

  std::vector<std::shared_ptr<arrow::Array>> arrays = {
      finish(idList),     finish(pdgList),     finish(massList),
      finish(energyList), finish(chargeList),  finish(vxList),
      finish(vyList),     finish(vzList),      finish(timeList),
      finish(pxList),     finish(pyList),      finish(pzList),
      finish(d0List),     finish(z0List),      finish(vprimList),
      finish(parentList), finish(primaryList),
  };

  auto table =
      arrow::Table::Make(ActsPlugins::ArrowUtil::particleSchema(), arrays);
  m_outputTable(ctx, ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
