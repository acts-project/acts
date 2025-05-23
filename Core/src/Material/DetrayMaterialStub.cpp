// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/DetrayExceptions.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

namespace Acts {

#define STUB_METHOD(type)                                       \
  std::unique_ptr<DetraySurfaceMaterial> type::toDetrayPayload( \
      const detray::io::volume_payload& /*volume*/) const {     \
    throw DetrayNotAvailableException();                        \
  }

// In STUB mode: all material related methods throw an exception to indicate
// that Detray is not available.

// clang-format off

STUB_METHOD(BinnedSurfaceMaterial)

template <>
STUB_METHOD(ProtoSurfaceMaterialT<Acts::BinUtility>)

template <>
STUB_METHOD(ProtoSurfaceMaterialT<std::vector<DirectedProtoAxis>>)

STUB_METHOD(detail::IGridSurfaceMaterialBase)

STUB_METHOD(HomogeneousSurfaceMaterial)

// clang-format on

}  // namespace Acts
