// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/GeoModel/detail/GeoShiftConverter.hpp"

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConversionError.hpp"
#include "ActsPlugins/GeoModel/detail/GeoBoxConverter.hpp"
#include "ActsPlugins/GeoModel/detail/GeoTrdConverter.hpp"
#include "ActsPlugins/GeoModel/detail/GeoTubeConverter.hpp"

#include <GeoModelHelpers/GeoShapeUtils.h>
#include <GeoModelKernel/GeoShapeShift.h>
#include <GeoModelKernel/GeoTube.h>

using namespace Acts;

namespace ActsPlugins::detail {

namespace {

template <typename ContainedShape, typename Converter>
Result<GeoModelSensitiveSurface> impl(const PVConstLink& geoPV,
                                      const GeoShapeShift& geoShift,
                                      const Transform3& absTransform,
                                      SurfaceBoundFactory& boundFactory,
                                      bool sensitive) {
  auto* trd = dynamic_cast<const ContainedShape*>(geoShift.getOp());

  if (trd == nullptr) {
    return GeoModelConversionError::WrongShapeForConverter;
  }

  return Converter{}(geoPV, *trd, absTransform * geoShift.getX(), boundFactory,
                     sensitive);
}

}  // namespace

Result<GeoModelSensitiveSurface> GeoShiftConverter::operator()(
    const PVConstLink& geoPV, const GeoShapeShift& geoShift,
    const Transform3& absTransform, SurfaceBoundFactory& boundFactory,
    bool sensitive) const {
  const auto opType = geoShift.getOp()->typeID();
  if (opType == GeoTrd::getClassTypeID()) {
    return impl<GeoTrd, detail::GeoTrdConverter>(geoPV, geoShift, absTransform,
                                                 boundFactory, sensitive);
  } else if (opType == GeoBox::getClassTypeID()) {
    return impl<GeoBox, detail::GeoBoxConverter>(geoPV, geoShift, absTransform,
                                                 boundFactory, sensitive);
  } else if (opType == GeoTube::getClassTypeID()) {
    return impl<GeoTube, detail::GeoTubeConverter>(
        geoPV, geoShift, absTransform, boundFactory, sensitive);
  }

  return GeoModelConversionError::WrongShapeForConverter;
}

}  // namespace ActsPlugins::detail
