/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/utils/vecmem_objects.hpp"

#include "get_queue.hpp"

// VecMem include(s).
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
#include <vecmem/memory/hip/device_memory_resource.hpp>
#include <vecmem/memory/hip/host_memory_resource.hpp>
#include <vecmem/memory/hip/managed_memory_resource.hpp>
#include <vecmem/utils/hip/async_copy.hpp>
#include <vecmem/utils/hip/copy.hpp>
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/memory/sycl/shared_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>
#include <vecmem/utils/sycl/copy.hpp>
#else
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#endif

namespace traccc::alpaka {

struct vecmem_objects::impl {

    /// Constructor
    explicit impl([[maybe_unused]] queue& q)
        :
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED) || defined(ALPAKA_ACC_GPU_HIP_ENABLED)
          m_host_mr(),
          m_device_mr(::alpaka::getNativeHandle(
              ::alpaka::getDev(details::get_queue(q)))),
          m_shared_mr(),
          m_copy(),
          m_async_copy(::alpaka::getNativeHandle(details::get_queue(q)))
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
          m_queue(::alpaka::getNativeHandle(details::get_queue(q))),
          m_host_mr(&m_queue),
          m_device_mr(&m_queue),
          m_shared_mr(&m_queue),
          m_copy(&m_queue),
          m_async_copy(&m_queue)
#else
          m_host_mr(),
          m_device_mr(),
          m_shared_mr(),
          m_copy(),
          m_async_copy()
#endif
    {
    }

#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    vecmem::cuda::host_memory_resource m_host_mr;
    vecmem::cuda::device_memory_resource m_device_mr;
    vecmem::cuda::managed_memory_resource m_shared_mr;
    vecmem::cuda::copy m_copy;
    vecmem::cuda::async_copy m_async_copy;
#elif defined(ALPAKA_ACC_GPU_HIP_ENABLED)
    vecmem::hip::host_memory_resource m_host_mr;
    vecmem::hip::device_memory_resource m_device_mr;
    vecmem::hip::managed_memory_resource m_shared_mr;
    vecmem::hip::copy m_copy;
    vecmem::hip::async_copy m_async_copy;
#elif defined(ALPAKA_ACC_SYCL_ENABLED)
    ::sycl::queue m_queue;
    vecmem::sycl::host_memory_resource m_host_mr;
    vecmem::sycl::device_memory_resource m_device_mr;
    vecmem::sycl::shared_memory_resource m_shared_mr;
    vecmem::sycl::copy m_copy;
    vecmem::sycl::async_copy m_async_copy;
#else
    vecmem::host_memory_resource m_host_mr;
    vecmem::host_memory_resource m_device_mr;
    vecmem::host_memory_resource m_shared_mr;
    vecmem::copy m_copy;
    vecmem::copy m_async_copy;
#endif

};  // struct vecmem_objects::impl

vecmem_objects::vecmem_objects(queue& q) : m_impl{std::make_unique<impl>(q)} {}

vecmem_objects::vecmem_objects(vecmem_objects&&) noexcept = default;

vecmem_objects::~vecmem_objects() = default;

vecmem_objects& vecmem_objects::operator=(vecmem_objects&&) noexcept = default;

vecmem::memory_resource& vecmem_objects::host_mr() const {
    return m_impl->m_host_mr;
}

vecmem::memory_resource& vecmem_objects::device_mr() const {
    return m_impl->m_device_mr;
}

vecmem::memory_resource& vecmem_objects::shared_mr() const {
    return m_impl->m_shared_mr;
}

vecmem::copy& vecmem_objects::copy() const {
    return m_impl->m_copy;
}

vecmem::copy& vecmem_objects::async_copy() const {
    return m_impl->m_async_copy;
}

}  // namespace traccc::alpaka
