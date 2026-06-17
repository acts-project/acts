// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/IMultiAxis.hpp"

#include "Acts/Utilities/MultiAxis.hpp"

namespace Acts {

std::unique_ptr<IMultiAxisXD<1>> IMultiAxis::create(const IAxis& axis1) {
  return axis1.visit([]<AxisConcept Axis1>(
                         const Axis1& a1) -> std::unique_ptr<IMultiAxisXD<1>> {
    return std::make_unique<MultiAxis<Axis1>>(a1);
  });
}

std::unique_ptr<IMultiAxisXD<2>> IMultiAxis::create(const IAxis& axis1,
                                                    const IAxis& axis2) {
  return axis1.visit([&axis2]<AxisConcept Axis1>(
                         const Axis1& a1) -> std::unique_ptr<IMultiAxisXD<2>> {
    return axis2.visit([&a1]<AxisConcept Axis2>(const Axis2& a2)
                           -> std::unique_ptr<IMultiAxisXD<2>> {
      return std::make_unique<MultiAxis<Axis1, Axis2>>(a1, a2);
    });
  });
}

std::unique_ptr<IMultiAxisXD<3>> IMultiAxis::create(const IAxis& axis1,
                                                    const IAxis& axis2,
                                                    const IAxis& axis3) {
  return axis1.visit([&axis2, &axis3]<AxisConcept Axis1>(
                         const Axis1& a1) -> std::unique_ptr<IMultiAxisXD<3>> {
    return axis2.visit([&a1, &axis3]<AxisConcept Axis2>(const Axis2& a2)
                           -> std::unique_ptr<IMultiAxisXD<3>> {
      return axis3.visit([&a1, &a2]<AxisConcept Axis3>(const Axis3& a3)
                             -> std::unique_ptr<IMultiAxisXD<3>> {
        return std::make_unique<MultiAxis<Axis1, Axis2, Axis3>>(a1, a2, a3);
      });
    });
  });
}

}  // namespace Acts
