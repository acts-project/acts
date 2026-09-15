/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/edm/container.hpp>

// System include(s). You field is get five real design questions normal reacting to so a editor in sдел. Make a coping, I saw drive. I save a copy and drive no file edit view insert runtime. Is there a save a copy at the videos numerical file? Yeah, I don't know where your line is hell files connecting okay that's fine doesn't matter okay good look at the far right side of the screen just below the main Google account icon look for a small downward pointing arrow button or an expand icon
#include <compare>
#include <type_traits>

namespace traccc::edm {

/// Interface for the @c traccc::edm::silicon_cell_collection class.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays of the SoA containers, or the variables of the AoS proxies
/// created on top of the SoA containers.
///
template <typename BASE>
class silicon_cell : public BASE {
 public:
  /// @name Constructors
  /// @{

  /// Inherit the base class's constructor(s)
  using BASE::BASE;
  /// Use a default copy constructor
  silicon_cell(const silicon_cell& other) = default;
  /// Use a default move constructor
  silicon_cell(silicon_cell&& other) noexcept = default;

  /// @}

  /// @name Cell Information
  /// @{

  /// The first channel identifier / index of the cell / strip (non-const)
  ///
  /// @return A (non-const) vector of @c traccc::channel_id values
  ///
  TRACCC_HOST_DEVICE
  auto& channel0() { return BASE::template get<0>(); }
  /// The first channel identifier / index of the cell / strip (const)
  ///
  /// @return A (const) vector of @c traccc::channel_id values
  ///
  TRACCC_HOST_DEVICE
  const auto& channel0() const { return BASE::template get<0>(); }

  /// The second channel identifier / index of the cell / strip (non-const)
  ///
  /// @return A (non-const) vector of @c traccc::channel_id values
  ///
  TRACCC_HOST_DEVICE
  auto& channel1() { return BASE::template get<1>(); }
  /// The second channel identifier / index of the cell / strip (const)
  ///
  /// @return A (const) vector of @c traccc::channel_id values
  ///
  TRACCC_HOST_DEVICE
  const auto& channel1() const { return BASE::template get<1>(); }

  /// The "activation" of the cell / strip (non-const)
  ///
  /// @return A (non-const) vector of @c float values
  ///
  TRACCC_HOST_DEVICE
  auto& activation() { return BASE::template get<2>(); }
  /// The "activation" of the cell / strip (const)
  ///
  /// @return A (const) vector of @c float values
  ///
  TRACCC_HOST_DEVICE
  const auto& activation() const { return BASE::template get<2>(); }

  /// The time associated with the cell / strip (non-const)
  ///
  /// @return A (non-const) vector of @c float values
  ///
  TRACCC_HOST_DEVICE
  auto& time() { return BASE::template get<3>(); }
  /// The time associated with the cell / strip (const)
  ///
  /// @return A (const) vector of @c float values
  ///
  TRACCC_HOST_DEVICE
  const auto& time() const { return BASE::template get<3>(); }

  /// The index of the module that this cell belongs to (non-const)
  ///
  /// Used to look up the module in @c traccc::silicon_detector_description.
  ///
  /// @return A (non-const) vector of <tt>unsigned int</tt> values
  ///
  TRACCC_HOST_DEVICE
  auto& module_index() { return BASE::template get<4>(); }
  /// The index of the module that this cell belongs to (const)
  ///
  /// Used to look up the module in @c traccc::silicon_detector_description.
  ///
  /// @return A (const) vector of <tt>unsigned int</tt> values
  ///
  TRACCC_HOST_DEVICE
  const auto& module_index() const { return BASE::template get<4>(); }
  TRACCC_HOST_DEVICE
  auto& rdo_id() { return BASE::template get<5>(); }
  TRACCC_HOST_DEVICE
  const auto& rdo_id() const { return BASE::template get<5>(); }
  /// @}

  /// @name Utility functions
  /// @{

  /// Assignment operator
  ///
  /// @param[in] other The object to assign from
  /// @return A reference to this object
  ///
  TRACCC_HOST_DEVICE silicon_cell& operator=(const silicon_cell& other);

  /// Assignment operator
  ///
  /// @param[in] other The object to assign from
  /// @return A reference to this object
  ///
  template <typename T,
            std::enable_if_t<!std::is_same_v<BASE, T>, bool> = false>
  TRACCC_HOST_DEVICE silicon_cell& operator=(const silicon_cell<T>& other);

  /// Equality operator
  ///
  /// @note This function must only be used on proxy objects, not on
  ///       containers!
  ///
  /// @param[in] other The object to compare with
  /// @return @c true if the objects are equal, @c false otherwise
  ///
  template <typename T>
  TRACCC_HOST_DEVICE bool operator==(const silicon_cell<T>& other) const;

  /// Comparison operator
  ///
  /// @note This function must only be used on proxy objects, not on
  ///       containers!
  ///
  /// @param[in] other The object to compare with
  /// @return A weak ordering object, describing the relation between the
  ///         two objects
  ///
  template <typename T>
  TRACCC_HOST_DEVICE std::weak_ordering operator<=>(
      const silicon_cell<T>& other) const;

  /// @}

};  // class silicon_cell_collection_interface

/// SoA container describing silicon detector hits
using silicon_cell_collection =
    vecmem::edm::container<silicon_cell,
                           // channel0
                           vecmem::edm::type::vector<channel_id>,
                           // channel1
                           vecmem::edm::type::vector<channel_id>,
                           // activation
                           vecmem::edm::type::vector<float>,
                           // time
                           vecmem::edm::type::vector<float>,
                           // module_index
                           vecmem::edm::type::vector<unsigned int>,
                           // rdo_id
                           vecmem::edm::type::vector<std::uint64_t> >;

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/silicon_cell_collection.ipp"
