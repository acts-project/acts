// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Material/MaterialSlab.hpp"

namespace Acts::detail {

class MaterialEffectsAccumulator {
 public:
  bool isVacuum() const { return m_accumulatedMaterial.isVacuum(); }

  void reset() { *this = MaterialEffectsAccumulator(); }

  void initialize(double maxXOverX0Step,
                  const ParticleHypothesis& particleHypothesis,
                  double initialMomentum);

  void accumulate(const Material& material, double pathLength, double qOverPin,
                  double qOverPout);

  std::optional<FreeMatrix> computeAdditionalFreeCovariance(
      const Vector3& direction);

 private:
  double m_maxXOverX0Step = 0;

  ParticleHypothesis m_particleHypothesis = ParticleHypothesis::pion();
  double m_initialMomentum = 0;

  MaterialSlab m_accumulatedMaterial;

  double m_varAngle = 0;
  double m_varPosition = 0;
  double m_covAnglePosition = 0;
};

}  // namespace Acts::detail
