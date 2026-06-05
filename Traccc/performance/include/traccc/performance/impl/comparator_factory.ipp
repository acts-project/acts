/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::details {

template <typename TYPE>
is_same_object<TYPE> comparator_factory<TYPE>::make_comparator(
    const TYPE& ref, scalar unc) const {

    return is_same_object<TYPE>{ref, unc};
}

}  // namespace traccc::details
