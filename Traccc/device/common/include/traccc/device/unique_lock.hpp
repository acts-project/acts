/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <mutex>

#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device {
/**
 * @brief An RAII-based lock, analogous to `std::unique_lock`.
 *
 * This unique lock type can be used to increase the safety of critical
 * sections at it guarantees that the lock will be released when it goes out of
 * scope.
 *
 * @tparam Mutex The underlying mutex type.
 *
 * @warning Be wary of using this code across different architectures, as
 * lockstep execution may seriously impact the behaviour of this code. In
 * particular, be wary of the fact that, on lockstep architectures, trying to
 * acquire this lock may result in deadlock as all locks will try to acquire
 * the lock at the same time. To resolve this, only try to acquire the lock
 * from one thread in a work group or block of threads. When using code written
 * in this way on non-lockstep architectures, be aware that the thread holding
 * the lock may exit a block before the other threads, releasing the lock early
 * and breaking thread safety. To avoid this problem, always synchronize
 * threads after the critical section.
 */
template <typename Mutex>
class unique_lock {
    public:
    using mutex_type = Mutex;

    /**
     * @brief Construct a unique lock without locking.
     */
    TRACCC_HOST_DEVICE
    unique_lock(mutex_type& m, std::defer_lock_t);

    /**
     * @brief Construct a unique lock, attempting to lock it.
     *
     * @warning This function returning does _not_ guarantee that the lock has
     * been acquired.
     *
     * @note Despite the warnings about acquiring locks in lockstep
     * architectures, calling this function across multiple threads is safe as
     * it is a non-blocking lock, e.g. it will fail for all but at most one
     * thread.
     */
    TRACCC_HOST_DEVICE
    unique_lock(mutex_type& m, std::try_to_lock_t);

    /**
     * @brief Construct a unique lock which was previously locked.
     */
    TRACCC_HOST_DEVICE
    unique_lock(mutex_type& m, std::adopt_lock_t);

    /**
     * @brief Destroy a lock, freeing the underlying mutex.
     */
    TRACCC_HOST_DEVICE
    ~unique_lock();

    /**
     * @brief Lock the lock, blocking until the operation succeeds.
     *
     * @warning On lockstep architectures, calling this method on more than a
     * single thread in a group will result in deadlock.
     */
    TRACCC_HOST_DEVICE
    void lock();

    /**
     * @brief Try to lock the lock without blocking.
     *
     * @note Calling this method from multiple threads in the same block is
     * safe.
     */
    TRACCC_HOST_DEVICE
    bool try_lock();

    /**
     * @brief Explicitly unlock the underlying lock.
     *
     * @warning Calling this method on a lock which has not been acquired
     * constitutes undefined behaviour.
     */
    TRACCC_HOST_DEVICE
    void unlock();

    /**
     * @brief Check if the lock is locked by this object.
     */
    TRACCC_HOST_DEVICE
    bool owns_lock() const;

    private:
    mutex_type* m_mutex_ptr = nullptr;
    bool m_owns_lock;
};
}  // namespace traccc::device

#include "impl/unique_lock.ipp"
