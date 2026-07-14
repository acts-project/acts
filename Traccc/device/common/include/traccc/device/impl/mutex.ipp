/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cassert>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {
template <typename T>
TRACCC_HOST_DEVICE mutex<T>::mutex(T& p) : m_atomic(p) {}

template <typename T>
TRACCC_HOST_DEVICE mutex<T>::mutex(const vecmem::device_atomic_ref<T>& r)
    : m_atomic(r) {}

template <typename T>
TRACCC_HOST_DEVICE void mutex<T>::lock() {
    while (!try_lock())
        ;
}

template <typename T>
TRACCC_HOST_DEVICE bool mutex<T>::try_lock() {
    assert(!m_is_locked);

    T false_v = static_cast<T>(false);
    bool s = m_atomic.compare_exchange_strong(false_v, static_cast<T>(true),
                                              vecmem::memory_order::acquire);

#ifndef NDEBUG
    m_is_locked |= s;
#endif

    return s;
}

template <typename T>
TRACCC_HOST_DEVICE void mutex<T>::unlock() {
    assert(m_is_locked);

    m_atomic.store(static_cast<T>(false), vecmem::memory_order::release);
}
}  // namespace traccc::device
