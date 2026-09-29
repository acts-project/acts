/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/global_index.hpp"
#include "traccc/edm/measurement_collection.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

/// Key type used for sorting the measurements.
///
/// A single 32 bit key, which is the measurement identifier assigned by the
/// clusterization as the index of the first cell of the cluster. The key is
/// unique for every measurement, and is sorted using a radix sort.
///
/// @note This relies on the input cells being ordered by surface identifier
///       so that ordering by cell index is the same as ordering by surface
///       identifier.
///
using measurement_sort_key_t = unsigned int;

/// Functor returning the sorting key of a measurement.
class measurement_sort_key_getter {
 public:
  /// Constructor with the (unsorted) measurements.
  ///
  /// @param measurements The view of the measurements to sort
  ///
  explicit TRACCC_HOST_DEVICE measurement_sort_key_getter(
      const edm::measurement_collection::const_view& measurements)
      : m_measurements(measurements) {}

  /// The operator returning one sorting key.
  ///
  /// @param index The index of the measurement
  /// @return The sorting key of the measurement
  ///
  TRACCC_HOST_DEVICE measurement_sort_key_t
  operator()(unsigned int index) const {
    const edm::measurement_collection::const_device measurements{
        m_measurements};
    return measurements.identifier().at(index);
  }

 private:
  /// The view of the (unsorted) measurements.
  edm::measurement_collection::const_view m_measurements;

};  // class measurement_sort_key_getter

/// Fill the sorting keys and the index sequence of the measurements.
///
/// @param[in]  globalIndex       The index of the current thread
/// @param[in]  measurements_view The unsorted measurements
/// @param[out] keys_view         The sorting keys of the measurements
/// @param[out] indices_view      The measurement indices to be sorted
///
TRACCC_HOST_DEVICE inline void fill_measurement_sort_keys(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& measurements_view,
    vecmem::data::vector_view<measurement_sort_key_t> keys_view,
    vecmem::data::vector_view<unsigned int> indices_view) {
  const edm::measurement_collection::const_device measurements{
      measurements_view};
  if (globalIndex >= measurements.size()) {
    return;
  }
  vecmem::device_vector<measurement_sort_key_t> keys{keys_view};
  vecmem::device_vector<unsigned int> indices{indices_view};
  keys.at(globalIndex) =
      measurement_sort_key_getter{measurements_view}(globalIndex);
  indices.at(globalIndex) = globalIndex;
}

/// Fill the output collection with the sorted measurements.
///
/// The identifier of every output measurement is set to its position in the
/// sorted collection. The cluster index is copied over unchanged, so that it
/// keeps pointing at the cluster that the measurement was created from.
///
/// @param[in]  globalIndex         The index of the current thread
/// @param[in]  input_view          The unsorted measurements
/// @param[out] output_view         The sorted measurements
/// @param[in]  sorted_indices_view The sorted measurement indices
///
TRACCC_HOST_DEVICE inline void fill_sorted_measurements(
    const global_index_t globalIndex,
    const edm::measurement_collection::const_view& input_view,
    edm::measurement_collection::view output_view,
    const vecmem::data::vector_view<const unsigned int>& sorted_indices_view) {
  const edm::measurement_collection::const_device input{input_view};
  if (globalIndex >= input.size()) {
    return;
  }
  edm::measurement_collection::device output{output_view};
  const vecmem::device_vector<const unsigned int> sorted_indices{
      sorted_indices_view};
  auto out = output.at(globalIndex);
  out = input.at(sorted_indices.at(globalIndex));
  out.identifier() = globalIndex;
}

}  // namespace traccc::device
