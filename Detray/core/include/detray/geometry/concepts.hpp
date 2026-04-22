// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/coordinates/coordinates.hpp"

// System include(s)
#include <concepts>
#include <ostream>
#include <type_traits>

namespace detray::concepts {

/// Concept for a type that performs coordinate shape transformations
template <class F>
concept coordinate_frame = requires() {
  requires concepts::point<typename F::loc_point>;
  requires concepts::algebra<typename F::algebra_type>;

  {
    F::global_to_local_3D(dtransform3D<typename F::algebra_type>(),
                          dpoint3D<typename F::algebra_type>(),
                          dvector3D<typename F::algebra_type>())
  } -> std::same_as<dpoint3D<typename F::algebra_type>>;

  {
    F::global_to_local(dtransform3D<typename F::algebra_type>(),
                       dpoint3D<typename F::algebra_type>(),
                       dvector3D<typename F::algebra_type>())
  } -> std::same_as<typename F::loc_point>;

  {
    F::local_to_global(dtransform3D<typename F::algebra_type>(),
                       dpoint3D<typename F::algebra_type>())
  } -> std::same_as<dpoint3D<typename F::algebra_type>>;

  {
    F::normal(dtransform3D<typename F::algebra_type>(),
              dpoint3D<typename F::algebra_type>())
  } -> std::same_as<dvector3D<typename F::algebra_type>>;
};

/// Concept for a geometrical shape
template <typename S, typename A>
concept shape = concepts::algebra<A> && requires(const S sh, std::ostream os) {
  requires std::is_enum_v<typename S::boundaries>;

  requires concepts::coordinate_frame<typename S::template local_frame_type<A>>;

  typename S::template result_type<dbool<A>>;
  typename S::template bounds_type<dscalar<A>>;

  S::dim;

  {
    sh.min_dist_to_boundary(typename S::template bounds_type<dscalar<A>>(),
                            dpoint2D<A>())
  } -> std::same_as<dscalar<A>>;

  {
    sh.template check_boundaries<A>(
        typename S::template bounds_type<dscalar<A>>(), dtransform3D<A>(),
        dpoint3D<A>(), dscalar<A>())
  } -> std::same_as<typename S::template result_type<dbool<A>>>;

  {
    sh.measure(typename S::template bounds_type<dscalar<A>>())
  } -> std::same_as<dscalar<A>>;

  {
    sh.merge(typename S::template bounds_type<dscalar<A>>(),
             typename S::template bounds_type<dscalar<A>>())
  } -> std::same_as<typename S::template bounds_type<dscalar<A>>>;

  {
    sh.template local_min_bounds<A>(
        typename S::template bounds_type<dscalar<A>>(), dscalar<A>())
  } -> std::same_as<darray<dscalar<A>, 6>>;

  {
    sh.template centroid<A>(typename S::template bounds_type<dscalar<A>>())
  } -> std::same_as<dpoint3D<A>>;

  {
    sh.template vertices<A>(typename S::template bounds_type<dscalar<A>>(),
                            dindex())
  } -> std::same_as<dvector<dpoint3D<A>>>;

  {
    sh.check_consistency(typename S::template bounds_type<dscalar<A>>(), os)
  } -> std::same_as<bool>;
};

/// Definition of a geometric object, that either has a shape or a local frame
/// @{
template <class G>
concept has_shape =
    concepts::algebra<typename G::algebra_type> &&
    concepts::shape<typename G::shape, typename G::algebra_type>;

template <class G>
concept has_local_frame = concepts::coordinate_frame<typename G::local_frame>;

template <class G>
concept geometric_object =
    concepts::algebra<typename G::algebra_type> &&
    (concepts::has_shape<G> || concepts::has_local_frame<G>);
/// @}

/// Frame types
/// @{
template <class F>
concept planar_frame =
    concepts::coordinate_frame<F> &&
    (std::same_as<F, cartesian2D<typename F::algebra_type>> ||
     std::same_as<F, polar2D<typename F::algebra_type>>);

template <class F>
concept rectilinear_frame =
    concepts::coordinate_frame<F> &&
    (std::same_as<F, cartesian2D<typename F::algebra_type>> ||
     std::same_as<F, cartesian3D<typename F::algebra_type>>);

template <class F>
concept cylindrical_frame =
    concepts::coordinate_frame<F> &&
    (std::same_as<F, cylindrical2D<typename F::algebra_type>> ||
     std::same_as<F, cylindrical3D<typename F::algebra_type>> ||
     std::same_as<F, concentric_cylindrical2D<typename F::algebra_type>>);

template <class F>
concept line_frame = concepts::coordinate_frame<F> &&
                     std::same_as<F, line2D<typename F::algebra_type>>;
/// @}

/// Shape types
/// @{
template <class S, class A>
concept planar_shape =
    concepts::shape<S, A> &&
    concepts::planar_frame<typename S::template local_frame_type<A>>;

template <class S, class A>
concept rectilinear_shape =
    concepts::shape<S, A> &&
    concepts::rectilinear_frame<typename S::template local_frame_type<A>>;

template <class S, class A>
concept cylindrical_shape =
    concepts::shape<S, A> &&
    concepts::cylindrical_frame<typename S::template local_frame_type<A>>;

template <class S, class A>
concept line_shape =
    concepts::shape<S, A> &&
    concepts::line_frame<typename S::template local_frame_type<A>>;
/// @}

/// Geometric object types
/// @{
template <class G>
concept planar_object =
    (concepts::geometric_object<G>) &&
    ((concepts::has_shape<G> &&
      concepts::planar_shape<typename G::shape, typename G::algebra_type>) ||
     (concepts::has_local_frame<G> &&
      concepts::planar_frame<typename G::local_frame>));

template <class G>
concept rectilinear_object =
    (concepts::geometric_object<G>) &&
    ((concepts::has_shape<G> &&
      concepts::rectilinear_shape<typename G::shape,
                                  typename G::algebra_type>) ||
     (concepts::has_local_frame<G> &&
      concepts::rectilinear_frame<typename G::local_frame>));

template <class G>
concept cylindrical_object =
    (concepts::geometric_object<G>) &&
    ((concepts::has_shape<G> &&
      concepts::cylindrical_shape<typename G::shape,
                                  typename G::algebra_type>) ||
     (concepts::has_local_frame<G> &&
      concepts::cylindrical_frame<typename G::local_frame>));

template <class G>
concept line_object =
    (concepts::geometric_object<G>) &&
    ((concepts::has_shape<G> &&
      concepts::line_shape<typename G::shape, typename G::algebra_type>) ||
     (concepts::has_local_frame<G> &&
      concepts::line_frame<typename G::local_frame>));
/// @}

template <class T>
concept planar = /*concepts::planar_shape<T, array<float>> ||*/
    concepts::planar_frame<T> || concepts::planar_object<T>;

template <class T>
concept rectilinear =
    /*concepts::rectilinear_shape<T, array<float>> ||*/
    concepts::rectilinear_frame<T> || concepts::rectilinear_object<T>;

template <class T>
concept cylindrical =
    /*concepts::cylindrical_shape<T, array<float>> ||*/
    concepts::cylindrical_frame<T> || concepts::cylindrical_object<T>;

template <class T>
concept line = /*concepts::line_shape<T, array<float>> ||*/
    concepts::line_frame<T> || concepts::line_object<T>;

}  // namespace detray::concepts
