// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"

#include <cmath>
#include <vector>

#include <boost/container/static_vector.hpp>

namespace ActsExamples {
using MuonSimHit = SimHit;
/// Container of space points.
using MuonSimHitContainer = std::vector<MuonSimHit>;
constexpr int g_fieldShift = 8;
// field translators
enum class MuonIdentifierFieldMaps {
  stationName = 40,
  stationEta = 32,
  stationPhi = 24,
  multilayer = 16,
  tubeLayer = 8,
  tube = 0,
};
struct muonMdtIdentifierFields {
  int8_t stationName = 0;
  int8_t stationEta = 0;
  int8_t stationPhi = 0;
  int8_t multilayer = 0;
  int8_t tubeLayer = 0;
  int8_t tube = 0;
};
muonMdtIdentifierFields splitId(Acts::GeometryIdentifier::Value theID) {
  muonMdtIdentifierFields f;
  f.tube = theID & 0xFF;
  theID = theID >> g_fieldShift;
  f.tubeLayer = theID & 0xFF;
  theID = theID >> g_fieldShift;
  f.multilayer = theID & 0xFF;
  theID = theID >> g_fieldShift;
  f.stationPhi = theID & 0xFF;
  theID = theID >> g_fieldShift;
  f.stationEta = theID & 0xFF;
  theID = theID >> g_fieldShift;
  f.stationName = theID & 0xFF;
  return f;
}
Acts::GeometryIdentifier::Value compressId(muonMdtIdentifierFields f) {
  Acts::GeometryIdentifier::Value out{0};
  out = out << g_fieldShift | f.stationName;
  out = out << g_fieldShift | f.stationEta;
  out = out << g_fieldShift | f.stationPhi;
  out = out << g_fieldShift | f.multilayer;
  out = out << g_fieldShift | f.tubeLayer;
  out = out << g_fieldShift | f.tube;
  return out;
}

}  // namespace ActsExamples
