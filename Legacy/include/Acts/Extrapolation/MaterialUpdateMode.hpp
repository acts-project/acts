// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
namespace Acts {

/// This is a steering enum to force the material update
/// it can be:
///  (1)  addNoise
/// (-1) removeNoise
/// Second is mainly for vertex reconstruction, but potentially dangeraous.

enum MaterialUpdateMode { addNoise = 1, removeNoise = -1 };
}