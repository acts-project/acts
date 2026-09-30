/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <fstream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace traccc::io::details {

/// Implementation detail for @c traccc::io::details::write_binary_soa
template <typename TYPE>
void write_binary_soa_variable(const TYPE& var, std::ostream& out_file) {
  // Make sure that the type works.
  static_assert(std::is_standard_layout_v<TYPE>,
                "Scalar type does not have a standard layout.");

  // Read the scalar variable as-is.
  out_file.write(reinterpret_cast<const char*>(&var), sizeof(TYPE));
}

/// Implementation detail for @c traccc::io::details::write_binary_soa
template <typename TYPE>
void write_binary_soa_variable(const vecmem::device_vector<TYPE>& var,
                               std::ostream& out_file) {
  // Make sure that the type works.
  static_assert(std::is_standard_layout_v<TYPE>,
                "Vector type does not have a standard layout.");

  // Write the size of the vector.
  const std::size_t size = var.size();
  out_file.write(reinterpret_cast<const char*>(&size), sizeof(std::size_t));

  // Write the contents of the vector.
  out_file.write(reinterpret_cast<const char*>(var.data()),
                 static_cast<std::streamsize>(size * sizeof(TYPE)));
}

/// Implementation detail for @c traccc::io::details::write_binary_soa
template <typename TYPE>
void write_binary_soa_variable(const vecmem::jagged_device_vector<TYPE>& var,
                               std::ostream& out_file) {
  // Make sure that the type works.
  static_assert(std::is_standard_layout_v<TYPE>,
                "Jagged vector type does not have a standard layout.");

  // Write the size of the "outer" vector.
  const std::size_t outer_size = var.size();
  out_file.write(reinterpret_cast<const char*>(&outer_size),
                 sizeof(std::size_t));

  // Write the sizes of the "inner" vectors.
  std::vector<std::size_t> inner_sizes(outer_size);
  for (std::size_t i = 0; i < outer_size; ++i) {
    inner_sizes.at(i) = var.at(i).size();
  }
  out_file.write(reinterpret_cast<const char*>(inner_sizes.data()),
                 outer_size * sizeof(typename std::size_t));

  // Write the inner vectors in multiple steps.
  for (std::size_t i = 0; i < outer_size; ++i) {
    out_file.write(reinterpret_cast<const char*>(var.at(i).data()),
                   inner_sizes.at(i) * sizeof(TYPE));
  }
}

/// Implementation detail for @c traccc::io::details::write_binary_soa
template <std::size_t INDEX, typename... VARTYPES,
          template <typename> class INTERFACE>
void write_binary_soa_impl(
    const vecmem::edm::device<vecmem::edm::schema<VARTYPES...>, INTERFACE>&
        container,
    std::ostream& out_file) {
  // Write the current variable.
  write_binary_soa_variable(container.template get<INDEX>(), out_file);

  // Recurse into the next variable.
  if constexpr (sizeof...(VARTYPES) > (INDEX + 1)) {
    write_binary_soa_impl<INDEX + 1>(container, out_file);
  }
}

/// Function reading an SoA container from a binary file
///
/// @param filename The full input filename
/// @param mr Is the memory resource to create the result container with
///
template <typename... VARTYPES, template <typename> class INTERFACE>
void write_binary_soa(
    std::string_view filename,
    const vecmem::edm::device<vecmem::edm::schema<VARTYPES...>, INTERFACE>&
        container) {
  // Open the output file.
  std::ofstream out_file(filename.data(), std::ios::binary);

  // Read all variables recursively.
  write_binary_soa_impl<0>(container, out_file);
}

}  // namespace traccc::io::details
