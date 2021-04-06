as// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once	

// CUDA plugin include(s)
#include "Acts/Plugins/Cuda/Utilities/NewMemoryManager/CudaStreamView.hpp"
#include "ErrorCheck.cuh"

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <atomic>
#include <cstddef>
#include <cstdint>

namespace Acts {
namespace Cuda {
namespace Nmm {

CudaStreamView::CudaStreamView(cudaStream_t stream) : stream_{stream} {}

cudaStream_t CudaStreamView::value() {
	return stream_;
}

CudaStreamView::cudaStream_t() {
	return CudaStreamView::value();
}

bool CudaStreamView::is_per_thread_default() {
	#ifdef CUDA_API_PER_THREAD_DEFAULT_STREAM
		return CudaStreamView::value() == cudaStreamPerThread || CudaStreamView::value() == 0;
	#else
		return CudaStreamView::value() == cudaStreamPerThread;
	#endif
}

bool CudaStreamView::is_default() {
	#ifdef CUDA_API_PER_THREAD_DEFAULT_STREAM
		return CudaStreamView::value() == cudaStreamPerLegacy;
	#else
		return CudaStreamView::value() == cudaStreamPerLegacy || CudaStreamView::value() == 0;
	#endif
}

void CudaStreamView::synchronize() {
	ACTS_CUDA_ERROR_CHECK(cudaStreamSynchronize(stream_));
}

/*
void CudaStreamView::synchronize_no_throw(){
	ACTS_CUDA_ERROR_CHECK(cudaStreamSynchronize(stream_));
}
*/

// Static CudaStreamView of the default stream (stream 0), for convenience
static constexpr CudaStreamView cuda_stream_default{};

// Static CudaStreamView of cudaStreamLegacy, for convenience
static CudaStreamView cuda_stream_legacy{cudaStreamLegacy};

// Static CudaStreamView of cudaStreamPerThread, for convenience
static CudaStreamView cuda_stream_per_thread{cudaStreamPerThread};

// Equality ciomparison operator for streams
// 
// @param[in] lhs the first stream view to compare
// @param[in] rhs the second stream view to compare
// @return true if equal, false if unequal
inline bool operator==(cuda_stream_view lhs, cuda_stream_view rhs) {
	return lhs.value() == rhs.value();
}

// Inequality comparison operator for streams
//
// @param[in] lhs the first stream view to compare
// @param[in] rhs the second stream view to compare
// @return true if unequal, false if equal
inline bool operator!=(cuda_stream_view lhs, cuda_stream_view rhs) { 
	return not(lhs == rhs); 
}

// Output stream operator for printing / logging streams
//
// @param[in] os the output ostream
// @param[in] sv the CudaStreamView to output
// @return std::ostream& the output ostream
inline std::ostream& operator<<(std::ostream& os, cuda_stream_view sv) {
	os << sv.value();
	return os;
}

} // namespace Nmm
} // namespace Cuda
} // namespace Acts