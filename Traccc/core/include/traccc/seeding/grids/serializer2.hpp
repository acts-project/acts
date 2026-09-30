/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc {

/** Struct that helps taking a two-dimensional binning into
 * a serial binning for data storage.
 *
 * Serializers allow to create a memory local environment if
 * advantageous.
 *
 **/
struct serializer2 {
  /** Create a serial bin from two individual bins
   *
   * @tparam faxis_t is the type of the first axis
   * @tparam saxis_t is the type of the second axis
   *
   * @param faxis the first axis
   * @param saxis the second axis, unused here
   * @param fbin first bin
   * @param sbin second bin
   *
   * @return a unsigned int for the memory storage
   */
  template <typename faxis_t, typename saxis_t>
  DETRAY_HOST_DEVICE unsigned int serialize(const faxis_t &faxis,
                                            const saxis_t & /*saxis*/,
                                            unsigned int fbin,
                                            unsigned int sbin) const {
    unsigned int offset = sbin * faxis.bins();
    return offset + fbin;
  }
};

}  // namespace traccc
