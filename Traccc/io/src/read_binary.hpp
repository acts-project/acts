/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// VecMem include(s).
#include <ios>
#include <vecmem/edm/host.hpp>
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstddef>
#include <fstream>
#include <string_view>
#include <type_traits>
#include <vector>

namespace traccc::io::details {

/// Function for reading a container from a binary file
///
/// @param filename The full input filename
/// @param mr Is the memory resource to create the result container with
///
/// TODO: Change container reading to not own its result object
template <typename container_t>
container_t read_binary_container(std::string_view filename,
                                  vecmem::memory_resource* mr = nullptr) {

    // Make sure that the chosen types work.
    static_assert(std::is_standard_layout_v<typename container_t::header_type>,
                  "Container header type must be standard layout.");
    static_assert(std::is_standard_layout_v<typename container_t::item_type>,
                  "Container item type must be standard layout.");

    // Open the input file.
    std::ifstream in_file(filename.data(), std::ios::binary);

    // Read the size of the header vector.
    std::size_t headers_size;
    in_file.read(reinterpret_cast<char*>(&headers_size), sizeof(std::size_t));

    // Read the sizes of the item vector.
    std::vector<std::size_t> items_size(headers_size);
    in_file.read(reinterpret_cast<char*>(items_size.data()),
                 static_cast<std::streamsize>(headers_size *
                                              sizeof(typename std::size_t)));

    // Create the result container, and set it to the correct (outer) size right
    // away.
    container_t result(headers_size, mr);

    // Read the header payload into memory.
    in_file.read(reinterpret_cast<char*>(result.get_headers().data()),
                 headers_size * sizeof(typename container_t::header_type));

    // Read the items in multiple steps.
    for (std::size_t i = 0; i < headers_size; ++i) {
        result.get_items().at(i).resize(items_size.at(i));
        in_file.read(
            reinterpret_cast<char*>(result.get_items().at(i).data()),
            items_size.at(i) * sizeof(typename container_t::item_type));
    }

    // Return the newly created container.
    return result;
}

/// Function for reading a collection from a binary file
///
/// @param filename The full input filename
/// @param mr Is the memory resource to create the result collection with
///
template <typename collection_t>
void read_binary_collection(collection_t& result, std::string_view filename) {

    // Make sure that the chosen types work.
    static_assert(std::is_standard_layout_v<typename collection_t::value_type>,
                  "Collection item type must be standard layout.");

    // Open the input file.
    std::ifstream in_file(filename.data(), std::ios::binary);

    // Read the size of the header vector.
    std::size_t size;
    in_file.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    // Set result to the correct size.
    result.resize(size);

    // Read the items into memory.
    in_file.read(reinterpret_cast<char*>(result.data()),
                 static_cast<std::streamsize>(
                     size * sizeof(typename collection_t::value_type)));
}

/// Implementation detail for @c traccc::io::details::read_binary_soa
template <typename TYPE>
void read_binary_soa_variable(TYPE& result, std::istream& in_file) {

    // Make sure that the type works.
    static_assert(std::is_standard_layout_v<TYPE>,
                  "Scalar type does not have a standard layout.");

    // Read the scalar variable as-is.
    in_file.read(reinterpret_cast<char*>(&result), sizeof(TYPE));
}

/// Implementation detail for @c traccc::io::details::read_binary_soa
template <typename TYPE, typename ALLOC>
void read_binary_soa_variable(std::vector<TYPE, ALLOC>& result,
                              std::istream& in_file) {

    // Make sure that the type works.
    static_assert(std::is_standard_layout_v<TYPE>,
                  "Vector type does not have a standard layout.");

    // Read the size of the vector.
    std::size_t size;
    in_file.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

    // Set result to the correct size.
    result.resize(size);

    // Read the items into memory.
    in_file.read(reinterpret_cast<char*>(result.data()),
                 static_cast<std::streamsize>(size * sizeof(TYPE)));
}

/// Implementation detail for @c traccc::io::details::read_binary_soa
template <typename TYPE, typename ALLOC1, typename ALLOC2>
void read_binary_soa_variable(
    std::vector<std::vector<TYPE, ALLOC2>, ALLOC1>& result,
    std::istream& in_file) {

    // Make sure that the type works.
    static_assert(std::is_standard_layout_v<TYPE>,
                  "Jagged vector type does not have a standard layout.");

    // Read the size of the "outer" vector.
    std::size_t outer_size;
    in_file.read(reinterpret_cast<char*>(&outer_size), sizeof(std::size_t));

    // Read the sizes of the "inner" vectors.
    std::vector<std::size_t> inner_sizes(outer_size);
    in_file.read(reinterpret_cast<char*>(inner_sizes.data()),
                 static_cast<std::streamsize>(outer_size *
                                              sizeof(typename std::size_t)));

    // Resize the outer vector.
    result.resize(outer_size);

    // Read the inner vectors in multiple steps.
    for (std::size_t i = 0; i < outer_size; ++i) {
        result.at(i).resize(inner_sizes.at(i));
        in_file.read(reinterpret_cast<char*>(result.at(i).data()),
                     inner_sizes.at(i) * sizeof(TYPE));
    }
}

/// Implementation detail for @c traccc::io::details::read_binary_soa
template <std::size_t INDEX, typename... VARTYPES,
          template <typename> class INTERFACE>
void read_binary_soa_impl(
    vecmem::edm::host<vecmem::edm::schema<VARTYPES...>, INTERFACE>& result,
    std::istream& in_file) {

    // Read the current variable.
    read_binary_soa_variable(result.template get<INDEX>(), in_file);

    // Recurse into the next variable.
    if constexpr (sizeof...(VARTYPES) > (INDEX + 1)) {
        read_binary_soa_impl<INDEX + 1>(result, in_file);
    }
}

/// Function reading an SoA container from a binary file
///
/// @param filename The full input filename
/// @param mr Is the memory resource to create the result container with
///
template <typename... VARTYPES, template <typename> class INTERFACE>
void read_binary_soa(
    vecmem::edm::host<vecmem::edm::schema<VARTYPES...>, INTERFACE>& result,
    std::string_view filename) {

    // Open the input file.
    std::ifstream in_file(filename.data(), std::ios::binary);

    // Read all variables recursively.
    read_binary_soa_impl<0>(result, in_file);
}

}  // namespace traccc::io::details
