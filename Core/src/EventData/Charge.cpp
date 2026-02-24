// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/Charge.hpp"

#include "Acts/EventData/ChargeConcept.hpp"
#include "Acts/Utilities/Diagnostics.hpp"

namespace Acts {

// ensure concrete classes satisfy the concepts

ACTS_PUSH_IGNORE_DEPRECATED()
static_assert(ChargeConcept<Neutral>, "Neutral does not fulfill ChargeConcept");
static_assert(ChargeConcept<SinglyCharged>,
              "SinglyCharged does not fulfill ChargeConcept");
static_assert(ChargeConcept<NonNeutralCharge>,
              "NonNeutralCharge does not fulfill ChargeConcept");
ACTS_POP_IGNORE_DEPRECATED()

static_assert(ChargeConcept<AnyCharge>,
              "AnyCharge does not fulfill ChargeConcept");

}  // namespace Acts
