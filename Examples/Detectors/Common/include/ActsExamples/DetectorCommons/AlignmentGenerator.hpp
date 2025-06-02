// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"

namespace ActsExamples::AlignmentGenerator {
    /// This generator does nothing 
    struct Nominal {
        /// @brief The call operation for the nominal alignment
        void operator()(Acts::Transform3& /*transform*/)  {}
            // No operation, this is a nominal alignment generator}
    };

    /// This generator applies a constant global shift to the transform
    struct GlobalShift {
        Acts::Vector3 shift;
       /// @brief The call operation applying the global shift
        void operator()(Acts::Transform3& transform) {
            transform.translation() += shift;
        }
    };

    /// This generator applies a constant global rotation of value `angle` around the `axis`
    struct GlobalRotation {
        /// The axis around which the rotation is applied
        Acts::Vector3 axis = Acts::Vector3::UnitZ();  ///< The rotation axis, default is Z-axis
        double angle = 0.0;  ///< The rotation angle in radians

        /// @brief The call operation applying the global rotation
        /// @param transform The transform to be rotated
        void operator()(Acts::Transform3& transform) {
            transform *= Acts::AngleAxis3(angle, axis);
        }
    };


}  // namespace ActsExamples