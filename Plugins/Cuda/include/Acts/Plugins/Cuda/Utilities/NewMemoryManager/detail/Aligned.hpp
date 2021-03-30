// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>

namespace Acts {
namespace Cuda {
namespace Nmm {
namespace detail {

//defines a default alignment used for a host memory allocated
static constexpr std::size_t NMM_DEFAULT_HOST_ALIGNMENT{alignof(std::max_align_t)};

//@return if `n` is a power of 2
constexpr bool is_pow2(std::size_t n) { return (0 == (n & (n - 1)));} 

//@return if `alignment` is a valid memory alignment
constexpr bool is_suppored_aligment(std::size_t alignment) { return is_pow2(alignment)}

//@param[in] v value to align
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return the aligned value, as one would except
constexpr std::size_t align_up(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(is_suppored_aligment(align_bytes));
	return (v + (align_bytes - 1)) & ~(align_bytes - 1);
}

//@param[in] v value to align
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return the aligned value, as one would except
constexpr std::size_t align_down(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(is_suppored_aligment(align_bytes));
	return v & ~(align_bytes - 1);
}

//@param[in] v value to check for alignment
//@param[in] alignment amount, in bytes, must be a power of 2
//
//@return true if aligned
constexpr bool is_aligned(std::size_t v, std::size_t align_bytes) noexcept {
	//if the alignment is not support, the program will end
	assert(is_suppored_aligment(align_bytes));
	return v == align_down(v, align_bytes);
}



} //namespace detail
} //namespace Nmm
} //namespace Acts
} //namespace Cuda