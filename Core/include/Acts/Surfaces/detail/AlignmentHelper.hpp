// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"

#include <tuple>

namespace Acts::detail {

// The container for derivative of local frame axis w.r.t. its
// rotation parameters. The first element is for x axis, second for y axis and
// last for z axis
using RotationToAxes =
    std::tuple<RotationMatrix3, RotationMatrix3, RotationMatrix3>;

/// @brief Evaluate the derivative of local frame axes vector w.r.t.
/// its rotation around local x/y/z axis
/// @Todo: add parameter for rotation axis order
///
/// @param compositeRotation The rotation that help places the composite object being rotated
/// @param relRotation The relative rotation of the surface with respect to the composite object being rotated
///
/// @return Derivative of local frame x/y/z axis vector w.r.t. its
/// rotation angles (extrinsic Euler angles) around local x/y/z axis
RotationToAxes rotationToLocalAxesDerivative(
    const RotationMatrix3& compositeRotation,
    const RotationMatrix3& relRotation = RotationMatrix3::Identity());

/// @brief Evaluate the Jacobian of the local-frame alignment parameters of a
/// component (e.g. a sensor) w.r.t. the local-frame alignment parameters of
/// a composite structure it belongs to (e.g. a stave or a layer).
///
/// Local-frame alignment parameters of an object with local-to-global
/// transform (R, c) are (dt, dw): a translation dt along the object's own
/// local axes, and small rotation angles dw about its own local axes, pivoting
/// on its local origin c. The moved object has R' = R * (1 + [dw]x) and
/// c' = c + R * dt, with [v]x the cross-product matrix of v.
///
/// A rigid motion (dT, dW) of the composite moves the component by
///
///     da/dA = | Rrel^T   -Rs^T [d]x Rc |
///             |   0         Rrel^T     |
///
/// with Rc, Rs the composite and component rotations, Rrel = Rc^T Rs the
/// component rotation relative to the composite, and d = cs - cc the
/// component origin relative to the composite origin in global coordinates.
/// The derivation is in docs/pages/alignment_composite_jacobians.md.
///
/// @param compositeTransform The local-to-global transform of the composite
/// @param componentTransform The local-to-global transform of the component
///
/// @return The 6x6 Jacobian d(dt, dw)_component / d(dT, dW)_composite
AlignmentMatrix compositeToComponentJacobian(
    const Transform3& compositeTransform, const Transform3& componentTransform);

/// @brief Evaluate the Jacobian of the ACTS surface alignment parameters
/// (see @c AlignmentIndices) w.r.t. the local-frame alignment parameters
/// defined in @c compositeToComponentJacobian, for the same surface.
///
/// The ACTS parameters are the translation of the surface center in global
/// coordinates and small rotations about the local axes, so the Jacobian is
/// diag(R, 1). This is the single place encoding that convention: multiplying
/// derivatives w.r.t. the ACTS alignment parameters by it gives derivatives
/// w.r.t. the local-frame parameters.
///
/// @param surfaceTransform The local-to-global transform of the surface
///
/// @return The 6x6 Jacobian d(ACTS parameters) / d(local-frame parameters)
AlignmentMatrix localFrameToAlignmentParametersJacobian(
    const Transform3& surfaceTransform);

}  // namespace Acts::detail
