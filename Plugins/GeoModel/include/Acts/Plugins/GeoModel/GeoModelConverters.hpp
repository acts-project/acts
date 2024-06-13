// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Plugins/GeoModel/GeoModelConversionError.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorElement.hpp"
#include "Acts/Plugins/GeoModel/IGeoShapeConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GenericGeoShapeConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoBoxConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoIntersectionAnnulusConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoShiftConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoTrdConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoTubeConverter.hpp"
#include "Acts/Plugins/GeoModel/detail/GeoUnionDoubleTrdConverter.hpp"
#include "Acts/Utilities/Result.hpp"

#include <memory>
#include <tuple>

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoLogVol.h>
#include <GeoModelKernel/GeoShape.h>

namespace Acts {

/// @brief The GeoBox converter
///
/// This is a dedicated converter for GeoBox shapes
using GeoBoxConverter =
    detail::GenericGeoShapeConverter<GeoBox, detail::GeoBoxConverter>;

/// @brief A dedicated converter for GeoInterseciton that describe annulus bounds
///
/// This is very much tailored to the AnnulusBounds class
using GeoIntersectionAnnulusConverter =
    detail::GenericGeoShapeConverter<GeoShapeIntersection,
                                     detail::GeoIntersectionAnnulusConverter>;

/// @brief The GeoShift + Trd/Box/Tube converter
///
/// This is a dedicated converter for GeoBox shapes
using GeoShiftConverter =
    detail::GenericGeoShapeConverter<GeoShapeShift, detail::GeoShiftConverter>;

/// @brief The GeoTrd converter
///
/// This is a dedicated converter for GeoTrd shapes
using GeoTrdConverter =
    detail::GenericGeoShapeConverter<GeoTrd, detail::GeoTrdConverter>;

/// @brief The GeoTube converter
///
/// This is a dedicated converter for GeoTube shapes
using GeoTubeConverter =
    detail::GenericGeoShapeConverter<GeoTube, detail::GeoTubeConverter>;

/// @brief The GeoTube converter
///
/// This is a dedicated converter for GeoTube shapes
using GeoUnionDoubleTrdConverter =
    detail::GenericGeoShapeConverter<GeoShapeUnion,
                                     detail::GeoUnionDoubleTrdConverter>;

}  // namespace Acts
