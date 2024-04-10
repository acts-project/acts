// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Material/MaterialValidater.hpp"

Acts::MaterialValidater::MaterialValidater(
    const Acts::MaterialValidater::Config& cfg,
    std::unique_ptr<const Acts::Logger> mlogger)
    : m_cfg(cfg), m_logger(std::move(mlogger)) {
  if (m_cfg.materialAssigner == nullptr) {
    throw std::invalid_argument("Missing material assigner");
  }
}

Acts::RecordedMaterialTrack Acts::MaterialValidater::recordMaterial(
    const GeometryContext& gctx, const MagneticFieldContext& mctx,
    const Vector3& position, const Vector3& direction) const {
  Acts::RecordedMaterialTrack mTrack;

  auto [surfaceAssignments, volumeAssignments] =
      m_cfg.materialAssigner->assignmentCandidates(gctx, mctx, position,
                                                   direction);

  return mTrack;
}