// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/IMultiAxis.hpp"

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/MultiAxis.hpp"

#include "detail/CreateMultiAxis3.hpp"

namespace Acts {

std::unique_ptr<IMultiAxis1D> IMultiAxis::create(const IAxis& axis1) {
  return axis1.visit(
      []<AxisConcept Axis1>(const Axis1& a1) -> std::unique_ptr<IMultiAxis1D> {
        return std::make_unique<MultiAxis<Axis1>>(a1);
      });
}

std::unique_ptr<IMultiAxis2D> IMultiAxis::create(const IAxis& axis1,
                                                 const IAxis& axis2) {
  return axis1.visit([&axis2]<AxisConcept Axis1>(
                         const Axis1& a1) -> std::unique_ptr<IMultiAxis2D> {
    return axis2.visit([&a1]<AxisConcept Axis2>(
                           const Axis2& a2) -> std::unique_ptr<IMultiAxis2D> {
      return std::make_unique<MultiAxis<Axis1, Axis2>>(a1, a2);
    });
  });
}

std::unique_ptr<IMultiAxis3D> IMultiAxis::create(const IAxis& axis1,
                                                 const IAxis& axis2,
                                                 const IAxis& axis3) {
  // The concrete first axis is resolved here; the (much larger) instantiation
  // of the remaining two axes lives in per-boundary-type translation units
  // (CreateMultiAxis3*.cpp) to keep this TU's compiler memory bounded. The
  // helper is only declared here, so no MultiAxis instantiation happens in
  // this TU.
  return axis1.visit([&axis2, &axis3]<AxisConcept Axis1>(
                         const Axis1& a1) -> std::unique_ptr<IMultiAxis3D> {
    return detail::createMultiAxis3<Axis1>(a1, axis2, axis3);
  });
}

}  // namespace Acts
