/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/container.hpp"

namespace traccc {

/// Header: unsigned int for the number of lin_circles per spacepoint bin

/// Item: transformed coordinate of doublet of middle-bottom or middle-top
struct lin_circle {
    // z origin
    scalar m_Zo;
    // cotangent of pitch angle
    scalar m_cotTheta;
    // reciprocal of square of distance between two spacepoints
    scalar m_iDeltaR;
    // error term for sp-pair without correlation of middle space point
    scalar m_Er;
    // u component in transformed coordinate
    scalar m_U;
    // v component in transformed coordinate
    scalar m_V;

    TRACCC_HOST_DEVICE
    const scalar& Zo() const { return m_Zo; }

    TRACCC_HOST_DEVICE
    const scalar& cotTheta() const { return m_cotTheta; }

    TRACCC_HOST_DEVICE
    const scalar& iDeltaR() const { return m_iDeltaR; }

    TRACCC_HOST_DEVICE
    const scalar& Er() const { return m_Er; }

    TRACCC_HOST_DEVICE
    const scalar& U() const { return m_U; }

    TRACCC_HOST_DEVICE
    const scalar& V() const { return m_V; }
};

/// Declare all lin_circle collection types
using lin_circle_collection_types = collection_types<lin_circle>;

}  // namespace traccc
