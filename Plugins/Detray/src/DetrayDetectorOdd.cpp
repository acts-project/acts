// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Per-metadata instantiation of the heavy detray detector code (geometry
// conversion and detector I/O) for the ODD metadata.

#include "ActsPlugins/Detray/DetrayDetectorIO.ipp"
#include "ActsPlugins/Detray/DetrayGeometryConverter.ipp"

namespace ActsPlugins {

ACTS_DETRAY_CONVERT_INSTANTIATION(DetrayMetadata::Odd);
ACTS_DETRAY_IO_INSTANTIATION(, DetrayMetadata::Odd)

}  // namespace ActsPlugins
