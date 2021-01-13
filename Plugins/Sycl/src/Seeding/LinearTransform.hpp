// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// SYCL plugin include(s).
#include "Acts/Plugins/Sycl/Seeding/detail/Types.hpp"

#include "../Utilities/Arrays.hpp"
#include "SpacePointType.hpp"

// SYCL include(s).
#include <CL/sycl.hpp>

// System include(s).
#include <cassert>
#include <cstdint>

namespace Acts::Sycl::detail {

/// Functor performing a linear coordinate transformation on spacepoint pairs
template <SpacePointType OtherSPType>
class LinearTransform {
  // Sanity check(s).
  static_assert((OtherSPType == SpacePointType::Bottom) ||
                    (OtherSPType == SpacePointType::Top),
                "Class must be instantiated with either "
                "Acts::Sycl::detail::SpacePointType::Bottom or "
                "Acts::Sycl::detail::SpacePointType::Top");

 public:
  /// Constructor with all the necessary arguments
  LinearTransform(uint32_t nMiddleSPs,
                  const device_array<DeviceSpacePoint>& middleSPs,
                  uint32_t nOtherSPs,
                  const device_array<DeviceSpacePoint>& otherSPs,
                  const device_array<uint32_t>& middleIndexLUT,
                  const device_array<uint32_t>& otherIndexLUT, uint32_t nEdges,
                  device_array<DeviceLinEqCircle>& resultArray)
      : m_nMiddleSPs(nMiddleSPs),
        m_middleSPs(middleSPs.get()),
        m_nOtherSPs(nOtherSPs),
        m_otherSPs(otherSPs.get()),
        m_middleIndexLUT(middleIndexLUT.get()),
        m_otherIndexLUT(otherIndexLUT.get()),
        m_nEdges(nEdges),
        m_resultArray(resultArray.get()) {}

  /// Operator performing the coordinate linear transformation
  void operator()(cl::sycl::nd_item<1> item) const {
    // Get the index to operate on.
    const auto idx = item.get_global_linear_id();
    if (idx >= m_nEdges) {
      return;
    }

    // Translate this one index into indices in the spacepoint arrays.
    // Note that using asserts with the CUDA backend of dpc++ is not working
    // quite correctly at the moment. :-( So these checks may need to be
    // disabled if you need to build for an NVidia backend in Debug mode.
    const uint32_t middleIndex = m_middleIndexLUT[idx];
    assert(middleIndex < m_nMiddleSPs);
    (void)m_nMiddleSPs;
    const uint32_t otherIndex = m_otherIndexLUT[idx];
    assert(otherIndex < m_nOtherSPs);
    (void)m_nOtherSPs;

    // Create a copy of the spacepoint objects for the current thread. On
    // dedicated GPUs this provides a better performance than accessing
    // variables one-by-one from global device memory.
    const DeviceSpacePoint middleSP = m_middleSPs[middleIndex];
    const DeviceSpacePoint otherSP = m_otherSPs[otherIndex];

    // Calculate some "helper variables" for the coordinate linear
    // transformation.
    const float cosPhiM = middleSP.x / middleSP.r;
    const float sinPhiM = middleSP.y / middleSP.r;

    const float deltaX = otherSP.x - middleSP.x;
    const float deltaY = otherSP.y - middleSP.y;
    const float deltaZ = otherSP.z - middleSP.z;

    const float x = deltaX * cosPhiM + deltaY * sinPhiM;
    const float y = deltaY * cosPhiM - deltaX * sinPhiM;
    const float iDeltaR2 = 1.f / (deltaX * deltaX + deltaY * deltaY);

    // Create the result object.
    DeviceLinEqCircle result;
    result.iDeltaR = cl::sycl::sqrt(iDeltaR2);
    result.cotTheta = deltaZ * result.iDeltaR;
    if constexpr (OtherSPType == SpacePointType::Bottom) {
      result.cotTheta = -(result.cotTheta);
    }
    result.zo = middleSP.z - middleSP.r * result.cotTheta;
    result.u = x * iDeltaR2;
    result.v = y * iDeltaR2;
    result.er =
        ((middleSP.varZ + otherSP.varZ) +
         (result.cotTheta * result.cotTheta) * (middleSP.varR + otherSP.varR)) *
        iDeltaR2;

    // Store the result object in device global memory.
    m_resultArray[idx] = result;
    return;
  }

 private:
  /// Total number of middle spacepoints
  uint32_t m_nMiddleSPs;
  /// Pointer to the middle spacepoints (in global device memory)
  const DeviceSpacePoint* m_middleSPs;
  /// Total number of "other" (bottom or top) spacepoints
  uint32_t m_nOtherSPs;
  /// Pointer to the "other" (bottom or top) spacepoints (in global device mem.)
  const DeviceSpacePoint* m_otherSPs;

  /// Look-Up Table from the iteration index to the middle spacepoint index
  const uint32_t* m_middleIndexLUT;
  /// Loop-Up Table from the iteration index to the "other" spacepoint index
  const uint32_t* m_otherIndexLUT;

  /// Total number of elements in the result array
  uint32_t m_nEdges;

  /// The result array in device global memory
  DeviceLinEqCircle* m_resultArray;

};  // class LinearTransform

}  // namespace Acts::Sycl::detail
