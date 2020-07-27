// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

/// Set up the CUDA device with the requested ID for the test
///
/// @param deviceID The integer ID of the device to use for the test
/// @param maxThreadsPerBlock The maximum number of threads supported by the
///        device per block.
void setupCudaDevice(int deviceID, int& maxThreadsPerBlock);
