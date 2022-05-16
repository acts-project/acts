// This file is part of the Acts project.
//
// Copyright (C) 2017-2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Hit.hpp"

#include "edm4hep/SimTrackerHit.h"

namespace Acts {

ActsFatras::Hit convertEDM4hepSimHit(const edm4hep::SimTrackerHit& sth);

}  // namespace Acts
