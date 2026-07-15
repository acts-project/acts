/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc::details {

template <typename T>
is_same_object<T>::is_same_object(const T& ref, scalar) : m_ref(ref) {}

template <typename T>
bool is_same_object<T>::operator()(const T& obj) const {
    return (obj == m_ref.get());
}

}  // namespace traccc::details
