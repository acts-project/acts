// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <iosfwd>

namespace ActsExamples::ParticleId {

bool isHadron(int pdg);
bool isLepton(int pdg);

enum class HadronType {
  Hadron = 1,
  BBbarMeson = 2,
  CCbarMeson = 3,
  BottomMeson = 4,
  BottomBaryon = 5,
  CharmedMeson = 6,
  CharmedBaryon = 7,
  StrangeMeson = 8,
  StrangeBaryon = 9,
  LightMeson = 10,
  LightBaryon = 11,
  Unknown = 12
};

std::ostream& operator<<(std::ostream& os, HadronType hadron);

HadronType hadronLabel(int pdg);

}  // namespace ActsExamples::ParticleId
