// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Axis.hpp"
#include "Acts/Utilities/MultiAxis.hpp"

#include "CreateMultiAxis3.hpp"

namespace Acts::detail {

template <typename FirstAxis>
std::unique_ptr<IMultiAxis3D> createMultiAxis3(const FirstAxis& first,
                                               const IAxis& second,
                                               const IAxis& third) {
  return second.visit([&first, &third]<AxisConcept Axis2>(
                          const Axis2& a2) -> std::unique_ptr<IMultiAxis3D> {
    return third.visit([&first, &a2]<AxisConcept Axis3>(
                           const Axis3& a3) -> std::unique_ptr<IMultiAxis3D> {
      return std::make_unique<MultiAxis<FirstAxis, Axis2, Axis3>>(first, a2,
                                                                  a3);
    });
  });
}

}  // namespace Acts::detail
