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

Line_t Seeder::constructLine(const double theta, const double y0,
                             Line_t::ParamVector patternPars) const {
  Line_t::ParamVector linePars{};
  linePars[toUnderlying(Line_t::ParIndex::y0)] = y0;
  linePars[toUnderlying(Line_t::ParIndex::x0)] =
      patternPars[toUnderlying(Line_t::ParIndex::x0)];

  Vector3 patternDir{makeDirectionFromPhiTheta(
      patternPars[toUnderlying(Line_t::ParIndex::phi)],
      patternPars[toUnderlying(Line_t::ParIndex::theta)])};

  double patternTanAlpha = patternDir.x() / patternDir.z();
  if (patternTanAlpha > std::numeric_limits<double>::epsilon()) {
    const Vector3 dir =
        makeDirectionFromAxisTangents(patternTanAlpha, tan(theta));
    linePars[toUnderlying(Line_t::ParIndex::phi)] =
        Acts::VectorHelpers::phi(dir);
    linePars[toUnderlying(Line_t::ParIndex::theta)] =
        Acts::VectorHelpers::theta(dir);
  } else {
    linePars[toUnderlying(Line_t::ParIndex::phi)] = 90._degree;
    linePars[toUnderlying(Line_t::ParIndex::theta)] = theta;
  }
  return Line_t(linePars);
}

std::ostream& CompositeSpacePointLineSeeder::SeedParameters::print(
    std::ostream& os) const {
  os << "SeedSolution: theta = " << theta / UnitConstants::degree
     << " deg, y0 = " << y0 << " +/- " << dY0
     << ", dTheta = " << dTheta / UnitConstants::mrad
     << " mrad, nStrawHits = " << nStrawHits;
  return os;
}

bool CompositeSpacePointLineSeeder::isValidLine(SeedParameters seedSol) const {
  if (seedSol.theta < m_cfg.thetaRange[0] ||
      seedSol.theta > m_cfg.thetaRange[1]) {
    return false;
  }
  if (seedSol.y0 < m_cfg.interceptRange[0] ||
      seedSol.y0 > m_cfg.interceptRange[1]) {
    return false;
  }
  return true;
}
