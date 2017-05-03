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
#include "ACTS/Concepts/AnyGrid.hpp"
#include "ACTS/Utilities/Definitions.hpp"

BOOST_TYPE_ERASURE_MEMBER((Acts)(concept)(detail)(has_getField), getField, 1);

namespace Acts {

namespace concept {

  namespace detail {

    namespace mpl = boost::mpl;

    using field_mapper_concept
        = mpl::vector<has_getField<Vector3D(const Vector3D&), const bte::_self>,
                      bte::copy_constructible<>,
                      bte::relaxed>;
  }  // namespace detail

  namespace bte = boost::type_erasure;

  template <typename T = bte::_self>
  using AnyFieldMapper = bte::any<detail::field_mapper_concept, T>;

}  // namespace concept
}  // namespace Acts
