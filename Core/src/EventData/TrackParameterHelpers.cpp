// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameterHelpers.hpp"

#include "Acts/Utilities/VectorHelpers.hpp"

bool Acts::isBoundVectorValid(const BoundVector& v, bool validateAngleRange,
                              double epsilon, double maxAbsEta) {
  return BoundParametrization::globalSpherical().isValid(v, validateAngleRange,
                                                         epsilon, maxAbsEta);
}

bool Acts::isFreeVectorValid(const FreeVector& v, double epsilon,
                             double maxAbsEta) {
  bool pos0Valid = std::isfinite(v[eFreePos0]);
  bool pos1Valid = std::isfinite(v[eFreePos1]);
  bool pos2Valid = std::isfinite(v[eFreePos2]);
  bool dir0Valid = std::isfinite(v[eFreeDir0]);
  bool dir1Valid = std::isfinite(v[eFreeDir1]);
  bool dir2Valid = std::isfinite(v[eFreeDir2]);
  bool dirValid = (std::abs(v.segment<3>(eFreeDir0).norm() - 1) - epsilon <= 0);
  bool qOverPValid = std::isfinite(v[eFreeQOverP]);
  bool timeValid = std::isfinite(v[eFreeTime]);

  bool etaValid = true;
  if (std::isfinite(maxAbsEta)) {
    double eta = VectorHelpers::eta(v.segment<3>(eFreeDir0));
    etaValid = std::isfinite(eta) && (std::abs(eta) - epsilon <= maxAbsEta);
  }

  return pos0Valid && pos1Valid && pos2Valid && dir0Valid && dir1Valid &&
         dir2Valid && dirValid && qOverPValid && timeValid && etaValid;
}
