// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameterHelpers.hpp"

bool Acts::isBoundVectorValid(const BoundVector& v, double epsilon) {
  auto pi = std::numbers::pi_v<double>;

  bool loc0Valid = std::isfinite(v[eBoundLoc0]);
  bool loc1Valid = std::isfinite(v[eBoundLoc1]);
  bool phiValid = std::isfinite(v[eBoundPhi]) &&
                  (v[eBoundPhi] + epsilon >= -pi) &&
                  (v[eBoundPhi] - epsilon < pi);
  bool thetaValid = std::isfinite(v[eBoundTheta]) &&
                    (v[eBoundTheta] + epsilon >= 0) &&
                    (v[eBoundTheta] - epsilon <= pi);
  bool qOverPValid = std::isfinite(v[eBoundQOverP]) &&
                     (std::abs(v[eBoundQOverP]) + epsilon > 0);
  bool timeValid = std::isfinite(v[eBoundTime]);

  return loc0Valid && loc1Valid && phiValid && thetaValid && qOverPValid &&
         timeValid;
}

bool Acts::isFreeVectorValid(const FreeVector& v, double epsilon) {
  bool pos0Valid = std::isfinite(v[eFreePos0]);
  bool pos1Valid = std::isfinite(v[eFreePos1]);
  bool pos2Valid = std::isfinite(v[eFreePos2]);
  bool dir0Valid = std::isfinite(v[eFreeDir0]);
  bool dir1Valid = std::isfinite(v[eFreeDir1]);
  bool dir2Valid = std::isfinite(v[eFreeDir2]);
  bool dirValid = (v.segment<3>(eFreeDir0).norm() + epsilon > 0);
  bool qOverPValid =
      std::isfinite(v[eFreeQOverP]) && (std::abs(v[eFreeQOverP]) + epsilon > 0);
  bool timeValid = std::isfinite(v[eFreeTime]);

  return pos0Valid && pos1Valid && pos2Valid && dir0Valid && dir1Valid &&
         dir2Valid && dirValid && qOverPValid && timeValid;
}
