/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cstdint>
#include <vecmem/memory/device_atomic_ref.hpp>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {
/*
 * A mutex object over type T.
 *
 * @warning This class assumes that the value is written to _only_ by mutex
 * class objects. Writing to the given pointers or atomic references in any
 * other way is undefined behaviour. Furthermore, it is assumed that the
 * initial value of the underlying pointer is false.
 *
 * @warning This is a spinlock. Do not use when more efficient implementations
 * are available.
 */
template <typename T = uint32_t>
class mutex {
    public:
    /*
     * Construct a mutex from a pointer.
     */
    TRACCC_HOST_DEVICE
    mutex(T &);

    /*
     * Construct a mutex from a vecmem atomic reference.
     */
    TRACCC_HOST_DEVICE
    mutex(const vecmem::device_atomic_ref<T> &);

    /*
     * Attempt to acquire a lock on the mutex. This method spins until a lock
     * is acquired.
     *
     * @warning On lockstep devices, only one thread per thread group (e.g.
     * warp) should call this function!
     */
    TRACCC_HOST_DEVICE
    void lock();

    /*
     * Try to acquire a lock on the mutex, returning whether the operation
     * succeeded or not. */
    TRACCC_HOST_DEVICE
    bool try_lock();

    /*
     * Unlock the mutex.
     *
     * @warning Using this method on a mutex that is not locked is undefined
     * behaviour.
     */
    TRACCC_HOST_DEVICE
    void unlock();

    private:
    const vecmem::device_atomic_ref<T> m_atomic;

#ifndef NDEBUG
    bool m_is_locked = false;
#endif
};
}  // namespace traccc::device

#include "impl/mutex.ipp"
