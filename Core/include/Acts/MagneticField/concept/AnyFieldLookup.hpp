// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/mpl/vector.hpp>
#include <boost/type_erasure/any.hpp>
#include <boost/type_erasure/builtin.hpp>
#include <boost/type_erasure/member.hpp>
#include <boost/type_erasure/relaxed.hpp>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getField), getField, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getFieldCell), getFieldCell, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_isInside), isInside, 1)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getNBins), getNBins, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getMin), getMin, 0)
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getMax), getMax, 0)
// clang-format on

namespace Acts {

namespace concept {

  /// @cond
  namespace bte = boost::type_erasure;

  namespace afl_detail {

    namespace mpl = boost::mpl;

    using field_cell_concept
        = mpl::vector<has_getField<Vector3D(const Vector3D&), const bte::_self>,
                      has_isInside<bool(const Vector3D&), const bte::_self>,
                      bte::copy_constructible<>,
                      bte::relaxed>;
  }  // namespace afl_detail
  /// @endcond

  /// @ingroup MagneticField
  /// @brief any-type for field cell interface
  ///
  /// @tparam T placeholder specifying how to store the underlying object
  ///
  /// @c any type for all copy-constructible types @c U providing the following
  /// interface:
  /// @code{.cpp}
  /// Acts::Vector3D U::getField(const Acts::Vector3D&) const;
  /// bool U::isInside(const Acts::Vector3D&) const;
  /// @endcode
  ///
  /// @note By default, the contained object is stored by-value (= copied) into
  /// the @c boost::type_erasure::any object. In order to store the value by (@c
  /// const) reference, pass <tt>(const) boost::type_erasure::_self&</tt> as
  /// template parameter.
  template <typename T = bte::_self>
  using AnyFieldCell   = bte::any<afl_detail::field_cell_concept, T>;

  /// @cond
  namespace afl_detail {
    using field_lookup_concept
        = mpl::vector<field_cell_concept,
                      has_getFieldCell<AnyFieldCell<>(const Vector3D&),
                                       const bte::_self>,
                      has_getNBins<std::vector<size_t>(), const bte::_self>,
                      has_getMin<std::vector<double>(), const bte::_self>,
                      has_getMax<std::vector<double>(), const bte::_self>>;
  }  // namespace afl_detail
  /// @endcond

  /// @ingroup MagneticField
  /// @brief any-type for field look-up interface
  ///
  /// @tparam T placeholder specifying how to store the underlying object
  ///
  /// This any-type is a specialization of the Acts::concept::AnyFieldCell
  /// concept. Every object of this type also implements the
  /// Acts::concept::AnyFieldCell concept.
  ///
  /// The exact interface required for wrapped objects of type @c U is the one
  /// specified in Acts::concept::AnyFieldCell augmented with the following
  /// methods:
  /// @code{.cpp}
  /// struct U {
  ///   // access the field cell at a given global position
  ///   Acts::concept::AnyFieldCell getFieldCell(const Acts::Vector3D&) const;
  ///
  ///  // access the number of bins of all axes of the grid
  ///   std::array<size_t, DIM> getNBins() const;
  ///
  ///  // access the minimum value of all axes of the grid
  ///   std::array<double, DIM> getMin() const;
  ///
  ///  // access the maximum value of all axes of the grid
  ///   std::array<double, DIM> getMax() const;
  /// }
  /// @endcode
  ///
  /// @note By default, the contained object is stored by-value (= copied) into
  /// the @c boost::type_erasure::any object. In order to store the value by (@c
  /// const) reference, pass <tt>(const) boost::type_erasure::_self&</tt> as
  /// template parameter.
  template <typename T = bte::_self>
  using AnyFieldLookup = bte::any<afl_detail::field_lookup_concept, T>;

}  // namespace concept
}  // namespace Acts
