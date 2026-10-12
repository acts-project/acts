// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

#include <edm4hep/CaloHitContribution.h>
#include <edm4hep/SimCalorimeterHit.h>

namespace ActsExamples {

/// One source per native simulated deposit, in the same collection order.
/// Retain the input frame and its collections when writing relations: handles
/// do not replace the persistent objects to which they refer.
struct EDM4hepCalorimeterSource {
  edm4hep::SimCalorimeterHit hit;
  edm4hep::CaloHitContribution contribution;
};

using EDM4hepCalorimeterSourceContainer = std::vector<EDM4hepCalorimeterSource>;

}  // namespace ActsExamples
