// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"

// Test include(s).
#include "detray/test/device/device_fixture.hpp"
#include "detray/test/framework/types.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/memory_resource.hpp>

// GoogleTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <cmath>
#include <memory>

namespace detray::test {

/// Data fixture for device test cases
template <detray::concepts::algebra A>
class transform_fixture : public device_fixture<dscalar<A>> {
  using scalar_t = dscalar<A>;
  using vector3_t = dvector3D<A>;

  using result_t = scalar_t;
  using base_fixture = device_fixture<result_t>;

 public:
  /// Constructor, setting up the inputs for all of the tests
  transform_fixture(vecmem::memory_resource& mr) : base_fixture(mr) {}

 protected:
  /// Function setting things up for a test
  virtual void SetUp() override {
    // Set up the result vectors
    base_fixture::SetUp();

    // Set up the input vectors.
    m_t1 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());
    m_t2 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());
    m_t3 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());

    m_v1 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());
    m_v2 = std::make_unique<vecmem::vector<vector3_t>>(this->size(),
                                                       &this->resource());

    // Initialise the input and output vectors.
    for (std::size_t i = 0; i < this->size(); ++i) {
      m_t1->at(i) = {static_cast<scalar_t>(1.1), static_cast<scalar_t>(2.2),
                     static_cast<scalar_t>(3.3)};
      m_t2->at(i) = {static_cast<scalar_t>(4.4), static_cast<scalar_t>(5.5),
                     static_cast<scalar_t>(6.6)};
      m_t3->at(i) = {static_cast<scalar_t>(7.7), static_cast<scalar_t>(8.8),
                     static_cast<scalar_t>(9.9)};

      m_v1->at(i) = {static_cast<scalar_t>(i * 0.6),
                     static_cast<scalar_t>((i + 1) * 1.2),
                     static_cast<scalar_t>((i + 2) * 1.3)};
      m_v2->at(i) = {static_cast<scalar_t>((i + 1) * 1.8),
                     static_cast<scalar_t>(i * 2.3),
                     static_cast<scalar_t>((i + 2) * 3.4)};
    }
  }

  /// Function tearing things down after the test
  virtual void TearDown() override {
    // Delete the vectors.
    m_t1.reset();
    m_t2.reset();
    m_t3.reset();
    m_v1.reset();
    m_v2.reset();

    // Tear down the base class.
    base_fixture::TearDown();
  }

  /// @name Inputs for the tests
  /// @{

  std::unique_ptr<vecmem::vector<vector3_t>> m_t1, m_t2, m_t3;
  std::unique_ptr<vecmem::vector<vector3_t>> m_v1, m_v2;

  /// @}

};  // class device_fixture

/// Functor running @c test_device_basics::transform3_ops
template <detray::concepts::algebra A>
class transform3_ops_functor {
  using scalar_t = dscalar<A>;
  using point3_t = dpoint3D<A>;
  using vector3_t = dvector3D<A>;
  using transform3_t = dtransform3D<A>;

 public:
  DETRAY_HOST_DEVICE void operator()(
      std::size_t i, vecmem::data::vector_view<const vector3_t> t1,
      vecmem::data::vector_view<const vector3_t> t2,
      vecmem::data::vector_view<const vector3_t> t3,
      vecmem::data::vector_view<const vector3_t> a,
      vecmem::data::vector_view<const vector3_t> b,
      vecmem::data::vector_view<scalar_t> output) const {
    // Create the VecMem vector(s).
    vecmem::device_vector<const vector3_t> vec_t1(t1), vec_t2(t2), vec_t3(t3),
        vec_a(a), vec_b(b);
    vecmem::device_vector<scalar_t> vec_output(output);

    // Perform the operation.
    auto ii = static_cast<typename decltype(vec_output)::size_type>(i);
    vec_output[ii] = transform3_ops(vec_t1[ii], vec_t2[ii], vec_t3[ii],
                                    vec_a[ii], vec_b[ii]);
  }

 private:
  /// Perform various operations using the @c transform3 type
  DETRAY_HOST_DEVICE
  scalar_t transform3_ops(vector3_t t1, vector3_t t2, vector3_t t3, vector3_t a,
                          vector3_t b) const {
    using namespace algebra;

    transform3_t tr1(t1, t2, t3);
    transform3_t tr2;
    tr2 = tr1;

    point3_t translation = tr2.translation();

    point3_t gpoint = tr2.point_to_global(a);
    point3_t lpoint = tr2.point_to_local(b);

    vector3_t gvec = tr2.vector_to_global(a);
    vector3_t lvec = tr2.vector_to_local(b);

    return {detray::vector::norm(translation) + detray::vector::perp(gpoint) +
            detray::vector::phi(lpoint) + detray::vector::dot(gvec, lvec)};
  }
};

}  // namespace detray::test
