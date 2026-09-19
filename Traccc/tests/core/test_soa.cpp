/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include <gtest/gtest.h>
#include <vecmem/memory/host_memory_resource.hpp>

#include "traccc/utils/soa.hpp"

using namespace traccc;

namespace {

using test_layout = soa_layout<std::uint32_t, float[3], std::uint32_t>;
using test_view = soa_view<test_layout>;
using test_buffer = soa_buffer<test_layout>;

/// A type that names the columns, in the shape the CKF wrappers will take.
struct named_array : private soa_device<test_layout> {
  using soa_device::soa_device;

  decltype(auto) first(unsigned int i) { return at<0>(i); }
  decltype(auto) first(unsigned int i) const { return at<0>(i); }
  decltype(auto) vec(unsigned int i, unsigned int j) { return at<1>(i, j); }
  decltype(auto) vec(unsigned int i, unsigned int j) const {
    return at<1>(i, j);
  }
  decltype(auto) last(unsigned int i) { return at<2>(i); }
  decltype(auto) last(unsigned int i) const { return at<2>(i); }
};

std::uintptr_t address(const void* p) {
  return reinterpret_cast<std::uintptr_t>(p);
}

}  // namespace

TEST(SoA, LayoutAlignment) {
  // Small columns are padded up to one sector.
  static_assert(test_layout::max_alignment == 32u);
  // An over-aligned column raises the alignment of every column.
  struct alignas(64) wide {
    float v[16];
  };
  static_assert(soa_layout<std::uint32_t, wide>::max_alignment == 64u);
}

TEST(SoA, LayoutOffsets) {
  static_assert(test_layout::total_size == 20u);
  static_assert(test_layout::offset_at<0> == 0u);
  static_assert(test_layout::offset_at<1> == 4u);
  static_assert(test_layout::offset_at<2> == 16u);
  static_assert(test_layout::size_at<1> == 12u);
  static_assert(test_layout::flat_size_at<1> == 4u);
  static_assert(test_layout::extent_at<1> == 3u);
  static_assert(std::is_same_v<test_layout::type_at<1>, float[3]>);
  static_assert(std::is_same_v<test_layout::flat_type_at<1>, float>);
}

TEST(SoA, LayoutPaddedSize) {
  static_assert(test_layout::padded_size(0u) == 0u);
  static_assert(test_layout::padded_size(1u) == 32u);
  static_assert(test_layout::padded_size(17u) == 32u);
  static_assert(test_layout::padded_size(32u) == 32u);
  static_assert(test_layout::padded_size(33u) == 64u);
}

TEST(SoA, ViewIsTrivial) {
  static_assert(std::is_trivially_copyable_v<test_view>);
  static_assert(std::is_trivially_default_constructible_v<test_view>);
  static_assert(sizeof(test_view) == 16u);

  vecmem::host_memory_resource mr;
  test_buffer b(3u, mr);

  const test_view v(b);
  ASSERT_EQ(v.ptr(), b.ptr());
  ASSERT_EQ(v.size(), 3u);
  ASSERT_EQ(v.capacity(), 32u);
}

TEST(SoA, BufferConstruct) {
  vecmem::host_memory_resource mr;
  test_buffer b(17u, mr);

  ASSERT_EQ(b.size(), 17u);
  ASSERT_NE(b.ptr(), nullptr);
  ASSERT_EQ(address(b.ptr()) % test_layout::max_alignment, 0u);
}

TEST(SoA, BufferEmpty) {
  vecmem::host_memory_resource mr;
  test_buffer b(0u, mr);
  ASSERT_EQ(b.size(), 0u);
  ASSERT_EQ(b.ptr(), nullptr);

  test_buffer d;
  ASSERT_EQ(d.size(), 0u);
  ASSERT_EQ(d.ptr(), nullptr);
}

TEST(SoA, BufferMove) {
  vecmem::host_memory_resource mr;
  test_buffer a(5u, mr);
  void* ptr = a.ptr();

  test_buffer b(std::move(a));
  ASSERT_EQ(b.ptr(), ptr);
  ASSERT_EQ(b.size(), 5u);

  test_buffer c;
  c = std::move(b);
  ASSERT_EQ(c.ptr(), ptr);
  ASSERT_EQ(c.size(), 5u);
}

TEST(SoA, BufferTooLarge) {
  vecmem::host_memory_resource mr;
  // 20 bytes per row: this many rows exceed the 32-bit byte range. The check
  // runs before the allocation, so no memory is requested.
  const unsigned int too_many =
      std::numeric_limits<unsigned int>::max() / 20u + 1u;
  ASSERT_THROW(test_buffer(too_many, mr), std::length_error);
  // A row count that wraps in padded_size is refused as well.
  ASSERT_THROW(test_buffer(std::numeric_limits<unsigned int>::max(), mr),
               std::length_error);
}

TEST(SoA, ScalarRoundTrip) {
  vecmem::host_memory_resource mr;
  test_buffer b(14u, mr);
  test_view v(b);
  soa_device<test_layout> d(v);
  const soa_device<test_layout>& cd = d;

  for (unsigned int i = 0u; i < 14u; ++i) {
    d.at<0>(i) = 2u * i;
    d.at<2>(i) = 1000u + i;
  }

  for (unsigned int i = 0u; i < 14u; ++i) {
    ASSERT_EQ(d.at<0>(i), 2u * i);
    ASSERT_EQ(d.at<2>(i), 1000u + i);
    ASSERT_EQ(cd.at<0>(i), 2u * i);
    ASSERT_EQ(cd.at<2>(i), 1000u + i);
  }

  // The const accessor returns by value, the other one by reference.
  static_assert(std::is_same_v<decltype(d.at<0>(0u)), std::uint32_t&>);
  static_assert(std::is_same_v<decltype(cd.at<0>(0u)), std::uint32_t>);
}

TEST(SoA, ArrayRoundTrip) {
  vecmem::host_memory_resource mr;
  test_buffer b(14u, mr);
  soa_device<test_layout> d(b);
  const soa_device<test_layout>& cd = d;

  // Write in reverse order, so that an overlap between columns shows up as a
  // clobbered value.
  for (unsigned int i = 14u; i-- > 0u;) {
    for (unsigned int j = 3u; j-- > 0u;) {
      d.at<1>(i, j) = static_cast<float>(100u * i + j) + 0.5f;
    }
  }

  for (unsigned int i = 0u; i < 14u; ++i) {
    for (unsigned int j = 0u; j < 3u; ++j) {
      ASSERT_EQ(d.at<1>(i, j), static_cast<float>(100u * i + j) + 0.5f);
      ASSERT_EQ(cd.at<1>(i, j), static_cast<float>(100u * i + j) + 0.5f);
    }
  }

  static_assert(std::is_same_v<decltype(d.at<1>(0u, 0u)), float&>);
  static_assert(std::is_same_v<decltype(cd.at<1>(0u, 0u)), float>);
}

TEST(SoA, CompileTimeElementIndex) {
  vecmem::host_memory_resource mr;
  test_buffer b(17u, mr);
  soa_device<test_layout> d(b);
  const soa_device<test_layout>& cd = d;

  for (unsigned int i = 0u; i < 17u; ++i) {
    d.at<1, 0>(i) = 100.f * static_cast<float>(i);
    d.at<1, 1>(i) = 100.f * static_cast<float>(i) + 1.f;
    d.at<1, 2>(i) = 100.f * static_cast<float>(i) + 2.f;
  }

  // The compile-time and the runtime element index must agree on the address
  // and on the value.
  for (unsigned int i = 0u; i < 17u; ++i) {
    ASSERT_EQ((&d.at<1, 0>(i)), &d.at<1>(i, 0u));
    ASSERT_EQ((&d.at<1, 1>(i)), &d.at<1>(i, 1u));
    ASSERT_EQ((&d.at<1, 2>(i)), &d.at<1>(i, 2u));
    ASSERT_EQ((cd.at<1, 1>(i)), cd.at<1>(i, 1u));
    ASSERT_EQ((cd.at<1, 2>(i)), 100.f * static_cast<float>(i) + 2.f);
  }

  static_assert(std::is_same_v<decltype(d.at<1, 0>(0u)), float&>);
  static_assert(std::is_same_v<decltype(cd.at<1, 0>(0u)), float>);
}

TEST(SoA, ColumnLayout) {
  vecmem::host_memory_resource mr;
  // 17 rows pad to 32, so the column stride is 32 elements.
  test_buffer b(17u, mr);
  soa_device<test_layout> d(b);
  const std::uintptr_t base = address(b.ptr());

  // Every column starts at a sector boundary.
  ASSERT_EQ((address(&d.at<0>(0u)) - base) % 32u, 0u);
  ASSERT_EQ((address(&d.at<1, 0>(0u)) - base) % 32u, 0u);
  ASSERT_EQ((address(&d.at<1, 1>(0u)) - base) % 32u, 0u);
  ASSERT_EQ((address(&d.at<1, 2>(0u)) - base) % 32u, 0u);
  ASSERT_EQ((address(&d.at<2>(0u)) - base) % 32u, 0u);

  // Columns follow each other at the padded stride, and a row is contiguous
  // within a column.
  ASSERT_EQ(address(&d.at<0>(0u)) - base, 0u);
  ASSERT_EQ(address(&d.at<1, 0>(0u)) - base, 4u * 32u);
  ASSERT_EQ(address(&d.at<1, 1>(0u)) - base, 8u * 32u);
  ASSERT_EQ(address(&d.at<1, 2>(0u)) - base, 12u * 32u);
  ASSERT_EQ(address(&d.at<2>(0u)) - base, 16u * 32u);
  ASSERT_EQ(address(&d.at<0>(1u)) - address(&d.at<0>(0u)), 4u);
  ASSERT_EQ(address(&d.at<1, 0>(16u)) - address(&d.at<1, 0>(0u)), 64u);
}

TEST(SoA, ColumnsAreDisjoint) {
  vecmem::host_memory_resource mr;
  test_buffer b(40u, mr);
  soa_device<test_layout> d(b);

  // Fill every column with its own pattern, then check that no column
  // overwrote another. 40 rows pad to 64, so the padding rows are in play.
  for (unsigned int i = 0u; i < 40u; ++i) {
    d.at<0>(i) = 1u;
    d.at<1, 0>(i) = 2.f;
    d.at<1, 1>(i) = 3.f;
    d.at<1, 2>(i) = 4.f;
    d.at<2>(i) = 5u;
  }

  for (unsigned int i = 0u; i < 40u; ++i) {
    ASSERT_EQ(d.at<0>(i), 1u);
    ASSERT_EQ((d.at<1, 0>(i)), 2.f);
    ASSERT_EQ((d.at<1, 1>(i)), 3.f);
    ASSERT_EQ((d.at<1, 2>(i)), 4.f);
    ASSERT_EQ(d.at<2>(i), 5u);
  }
}

TEST(SoA, NamedWrapper) {
  vecmem::host_memory_resource mr;
  test_buffer b(14u, mr);
  named_array arr(b);
  const named_array& carr = arr;

  // The wrapper adds nothing to the accessor.
  static_assert(sizeof(named_array) == sizeof(soa_device<test_layout>));

  for (unsigned int i = 0u; i < 14u; ++i) {
    arr.first(i) = 2u * i;
    arr.last(i) = 3u * i;
    for (unsigned int j = 0u; j < 3u; ++j) {
      arr.vec(i, j) = static_cast<float>(100u * i + j);
    }
  }

  for (unsigned int i = 0u; i < 14u; ++i) {
    ASSERT_EQ(carr.first(i), 2u * i);
    ASSERT_EQ(carr.last(i), 3u * i);
    for (unsigned int j = 0u; j < 3u; ++j) {
      ASSERT_EQ(carr.vec(i, j), static_cast<float>(100u * i + j));
    }
  }
}

TEST(SoA, TwoBuffersCopy) {
  vecmem::host_memory_resource mr;
  test_buffer in(20u, mr);
  test_buffer out(20u, mr);
  soa_device<test_layout> din(in);
  soa_device<test_layout> doubt(out);
  const soa_device<test_layout>& cin = din;

  for (unsigned int i = 0u; i < 20u; ++i) {
    din.at<0>(i) = i;
    din.at<1, 0>(i) = static_cast<float>(i);
    din.at<1, 1>(i) = static_cast<float>(i) + 0.25f;
    din.at<1, 2>(i) = static_cast<float>(i) + 0.5f;
    din.at<2>(i) = 100u + i;
  }

  // A row copy between two buffers, column by column, as condense_tracks
  // will do it.
  for (unsigned int i = 0u; i < 20u; ++i) {
    doubt.at<0>(i) = cin.at<0>(i);
    doubt.at<1, 0>(i) = cin.at<1, 0>(i);
    doubt.at<1, 1>(i) = cin.at<1, 1>(i);
    doubt.at<1, 2>(i) = cin.at<1, 2>(i);
    doubt.at<2>(i) = cin.at<2>(i);
  }

  for (unsigned int i = 0u; i < 20u; ++i) {
    ASSERT_EQ(doubt.at<0>(i), i);
    ASSERT_EQ((doubt.at<1, 0>(i)), static_cast<float>(i));
    ASSERT_EQ((doubt.at<1, 1>(i)), static_cast<float>(i) + 0.25f);
    ASSERT_EQ((doubt.at<1, 2>(i)), static_cast<float>(i) + 0.5f);
    ASSERT_EQ(doubt.at<2>(i), 100u + i);
  }
}

TEST(SoA, NonGlobalFlagOnHost) {
  // The flag only changes device code; on the host both instantiations
  // behave the same.
  vecmem::host_memory_resource mr;
  test_buffer b(4u, mr);
  soa_device<test_layout, false> d(b);
  d.at<0>(3u) = 7u;
  ASSERT_EQ(d.at<0>(3u), 7u);
}
