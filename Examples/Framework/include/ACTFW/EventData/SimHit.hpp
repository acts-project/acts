// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ActsFatras/EventData/Hit.hpp"

namespace FW {

using SimHit = ::ActsFatras::Hit;
/// Store hits ordered by geometry identifier.
using SimHitContainer = GeometryIdMultiset<::ActsFatras::Hit>;

}  // end of namespace FW
