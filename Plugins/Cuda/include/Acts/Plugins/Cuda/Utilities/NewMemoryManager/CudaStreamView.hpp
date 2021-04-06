// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once	

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <atomic>
#include <cstddef>
#include <cstdint>

namespace Acts {
namespace Cuda {
namespace Nmm {

	class CudaStreamView {
		public:
			constexpr CudaStreamView()                        = default;
			constexpr CudaStreamView(CudaStreamView const&) = default;
			constexpr CudaStreamView(CudaStreamView&&)      = default;
			constexpr CudaStreamView& operator=(CudaStreamView const&) = default;
			constexpr CudaStreamView& operator=(CudaStreamView&&) = default;
			~CudaStreamView()                                       = default;

			// Implicit conversion from cudaStream_t
			constexpr CudaStreamView(cudaStream_t stream);

			// Returns the wrppped stream
			constexpr cudaStream_t value() const noexcept;
	
			// Explicit conversion to cudaStream_t
			explicit constexpr operator cudaStream_t() const noexcept;

			// Return true if the wrapped stream is the CUDA per-thread default stream
			bool is_per_thread_default() const noexcept;

			// Return true if the wrapped stream is explicitly the CUDA legacy default stream
			bool is_default() const noexcept;

			// Synchronize the viewed CUDA stream
			void synchronize() const;

			// Synchronize the viewed CUDA stream, don't throw if there is an error
			void synchronize_no_throw() const noexcept;

		private:
			cudaStream_t stream_{0};
	};

} // namespace Nmm
} // namespace Cuda
} // namespace Acts