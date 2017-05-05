// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
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
#include "ACTS/Utilities/Definitions.hpp"

// clang-format off
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_getField), getField, 1);
BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(afl_detail)(has_isInside), isInside, 1);
// clang-format on

namespace Acts {

namespace concept {

  /// @cond
  namespace bte = boost::type_erasure;

  namespace afl_detail {

    namespace mpl = boost::mpl;

    using field_lookup_concept
        = mpl::vector<has_getField<Vector3D(const Vector3D&), const bte::_self>,
                      has_isInside<bool(const Vector3D&), const bte::_self>,
                      bte::copy_constructible<>,
                      bte::relaxed>;
  }  // namespace afl_detail
  /// @endcond

  /// @ingroup MagneticField
  /// @brief any-type for field look-up interface
  ///
  /// @tparam T placeholder specifying how the value is stored
  ///
  /// @c any type for all types @c U providing the following interface:
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
  using AnyFieldLookup = bte::any<afl_detail::field_lookup_concept, T>;

}  // namespace concept
}  // namespace Acts
