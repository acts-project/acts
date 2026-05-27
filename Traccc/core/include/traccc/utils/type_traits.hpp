/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::details {

/// Helper trait for detecting when a type is a non-const version of another
///
/// This comes into play multiple times to enable certain constructors
/// conditionally through SFINAE.
///
template <typename CTYPE, typename NCTYPE>
struct is_same_nc {
    static constexpr bool value = false;
};

template <typename TYPE>
struct is_same_nc<TYPE, TYPE> {
    static constexpr bool value = true;
};

template <typename TYPE>
struct is_same_nc<const TYPE, TYPE> {
    static constexpr bool value = true;
};

}  // namespace traccc::details
