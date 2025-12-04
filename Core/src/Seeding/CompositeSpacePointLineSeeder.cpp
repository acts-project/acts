// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"

using namespace Acts::Experimental;

using Seeder = CompositeSpacePointLineSeeder;
using namespace Acts::UnitLiterals;

Seeder::CompositeSpacePointLineSeeder(const Config& cfg,
                                      std::unique_ptr<const Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

CompositeSpacePointLineSeeder::SeedParam_t Seeder::constructLine(
    const double theta, const double y0, SeedParam_t patternParams) const {
  SeedParam_t lineParams{};
  lineParams[toUnderlying(Line_t::ParIndex::y0)] = y0;
  lineParams[toUnderlying(Line_t::ParIndex::x0)] =
      patternParams[toUnderlying(Line_t::ParIndex::x0)];
  Vector3 patternDir{makeDirectionFromPhiTheta(
      patternParams[toUnderlying(Line_t::ParIndex::phi)],
      patternParams[toUnderlying(Line_t::ParIndex::theta)])};

  double patternTanAlpha = patternDir.x() / patternDir.z();
  if (patternTanAlpha > std::numeric_limits<double>::epsilon()) {
    const Vector3 dir =
        makeDirectionFromAxisTangents(patternTanAlpha, tan(theta));
    lineParams[toUnderlying(Line_t::ParIndex::phi)] =
        Acts::VectorHelpers::phi(dir);
    lineParams[toUnderlying(Line_t::ParIndex::theta)] =
        Acts::VectorHelpers::theta(dir);
  } else {
    lineParams[toUnderlying(Line_t::ParIndex::phi)] = 90._degree;
    lineParams[toUnderlying(Line_t::ParIndex::theta)] = theta;
  }
  return lineParams;
}

std::ostream& CompositeSpacePointLineSeeder::SeedParameters::print(
    std::ostream& os) const {
  os << "SeedSolution: theta = " << theta / UnitConstants::degree
     << " deg, y0 = " << y0 << " +/- " << dY0
     << ", dTheta = " << dTheta / UnitConstants::mrad
     << " mrad, nStrawHits = " << nStrawHits;
  return os;
}

bool CompositeSpacePointLineSeeder::isValidLine(
    const SeedParameters& seedSol) const {
  if (m_cfg.noCutsOnSeedParams) {
    return true;
  }
  if (seedSol.theta < m_cfg.thetaRange[0] ||
      seedSol.theta > m_cfg.thetaRange[1]) {
    ACTS_DEBUG("rejecting theta " << seedSol.theta);
    return false;
  }
  if (seedSol.y0 < m_cfg.interceptRange[0] ||
      seedSol.y0 > m_cfg.interceptRange[1]) {
    ACTS_DEBUG("rejecting y0 " << seedSol.y0);
    return false;
  }
  return true;
}
