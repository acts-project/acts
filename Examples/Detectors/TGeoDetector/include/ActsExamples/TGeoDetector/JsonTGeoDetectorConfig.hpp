// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/TGeo/TGeoCylinderDiscSplitter.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

namespace ActsExamples {

/// Read/Write config for layer configuration triplets
template<typename T>
void from_json(const nlohmann::json& j,
               ActsExamples::TGeoDetector::Config::LayerTriplet<T>& ltr) {

      std::array<T, 3> valueTriplet;
      //valueTriplet = j.get<T>();
      if constexpr (std::is_same_v<T, bool>) {
        valueTriplet = j.get<std::array<bool, 3>>();
      }
      else if constexpr (std::is_same_v<T, double>) {
        valueTriplet = j.get<std::array<double, 3>>();
      }
      else if constexpr (std::is_same_v<T, std::string>) {
        valueTriplet = j.get<std::array<std::string, 3>>();
      }
      else if constexpr (std::is_same_v<T, std::vector<std::string>>) {
        valueTriplet = j.get<std::array<std::vector<std::string>, 3>>();
      }
      else if constexpr (std::is_same_v<T, std::pair<double, double>>) {
        valueTriplet = j.get<std::array<std::pair<double, double>, 3>>();
      }
      ltr.negative = valueTriplet[0];
      ltr.central  = valueTriplet[1];
      ltr.positive = valueTriplet[2];
}

template<typename T>
ActsExamples::TGeoDetector::Config::LayerTriplet<T> to_json(const nlohmann::json& j);

/// Read config for cylinder/disc module splitter
void from_json(const nlohmann::json& j,
               Acts::TGeoCylinderDiscSplitter::Config& msc) {
  /// Number of segments in phi for a disc
  msc.cylinderPhiSegments = j["geo-tgeo-cyl-nphi-segs"];
  /// Number of segments in r for a disk
  msc.cylinderLongitudinalSegments = j["geo-tgeo-cyl-nz-segs"];
  /// Number of segments in phi for a disc
  msc.discPhiSegments = j["geo-tgeo-disc-nphi-segs"];
  /// Number of segments in r for a disk
  msc.discRadialSegments = j["geo-tgeo-disc-nr-segs"];
}

Acts::TGeoCylinderDiscSplitter::Config to_json(const nlohmann::json& j);

}  // namespace ActsExamples