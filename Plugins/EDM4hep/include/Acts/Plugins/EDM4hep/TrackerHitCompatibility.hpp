// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Compatibility with EDM4hep < 0.99 and >= 0.99
#if __has_include(<edm4hep/TrackerHit3D.h>)
#include "edm4hep/TrackerHit3D.h"
#include "edm4hep/TrackerHit3DCollection.h"
#else
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
namespace edm4hep {
using TrackerHit3DCollection = edm4hep::TrackerHitCollection;
using TrackerHit3D = edm4hep::TrackerHit;
}  // namespace edm4hep
#endif
