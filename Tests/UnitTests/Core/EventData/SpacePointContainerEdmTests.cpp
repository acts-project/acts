// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Holders.hpp"

#include <any>
#include <cmath>
#include <iterator>
#include <ranges>
#include <vector>

using namespace Acts;

namespace ActsTests {

struct SpacePoint {
  SpacePoint() = default;
  SpacePoint(float ix, float iy, float iz, float ivarR, float ivarZ)
      : x(ix), y(iy), z(iz), varR(ivarR), varZ(ivarZ) {}
  float x{};
  float y{};
  float z{};
  float varR{};
  float varZ{};
};
using SpacePointCollection = std::vector<SpacePoint>;

class Adapter {
 public:
  friend SpacePointContainer<Adapter, detail::RefHolder>;
  friend SpacePointContainer<Adapter, detail::ValueHolder>;
  using value_type = SpacePoint;
  using ValueType = value_type;

  Adapter(SpacePointCollection&&) = delete;
  explicit Adapter(const SpacePointCollection& externalCollection)
      : m_storage(&externalCollection) {}

 private:
  std::size_t size_impl() const { return storage().size(); }

  float x_impl(std::size_t idx) const { return storage()[idx].x; };
  float y_impl(std::size_t idx) const { return storage()[idx].y; };
  float z_impl(std::size_t idx) const { return storage()[idx].z; };
  float varianceR_impl(std::size_t idx) const { return storage()[idx].varR; }
  float varianceZ_impl(std::size_t idx) const { return storage()[idx].varZ; }

  const SpacePoint& get_impl(std::size_t idx) const { return storage()[idx]; }

  std::any component_impl(HashedString key, std::size_t /*n*/) const {
    using namespace HashedStringLiteral;
    switch (key) {
      case "TopStripVector"_hash:
      case "BottomStripVector"_hash:
      case "StripCenterDistance"_hash:
      case "TopStripCenterPosition"_hash:
        return Vector3(0., 0., 0.);
      default:
        throw std::runtime_error("no such component " + std::to_string(key));
    }
  }

  const SpacePointCollection& storage() const { return *m_storage; }

 private:
  const SpacePointCollection* m_storage{};
};

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(spacepoint_container_edm_traits) {
  using adapter_t = Adapter;
  using container_t = SpacePointContainer<adapter_t, detail::RefHolder>;
  using proxy_t = SpacePointProxy<container_t>;
  using iterator_t =
      detail::ContainerIterator<container_t, proxy_t&, std::size_t, false>;

  static_assert(std::ranges::range<container_t>);
  static_assert(std::same_as<typename iterator_t::iterator_category,
                             std::random_access_iterator_tag>);
  static_assert(
      std::same_as<typename std::iterator_traits<iterator_t>::iterator_category,
                   std::random_access_iterator_tag>);
}

BOOST_AUTO_TEST_CASE(spacepoint_container_edm_constructors) {
  std::size_t nExternalPoints = 10;
  SpacePointCollection externalCollection(nExternalPoints);

  SpacePointContainerConfig spConfig;
  SpacePointContainerOptions spOptions;

  Adapter adapterForRef(externalCollection);
  SpacePointContainer<Adapter, detail::RefHolder> spContainerRef(
      spConfig, spOptions, adapterForRef);

  SpacePointContainer<Adapter, detail::ValueHolder> spContainerVal(
      spConfig, spOptions, Adapter(externalCollection));

  BOOST_CHECK_EQUAL(spContainerRef.size(), nExternalPoints);
  BOOST_CHECK_EQUAL(spContainerVal.size(), nExternalPoints);
}

BOOST_AUTO_TEST_CASE(spacepoint_container_edm_functionalities) {
  std::size_t nExternalPoints = 100;
  SpacePointCollection externalCollection;
  externalCollection.reserve(nExternalPoints);
  for (std::size_t i = 0; i < nExternalPoints; ++i) {
    externalCollection.emplace_back(1.f * i, 1.5f * i, 2.f * i, 2.5f * i,
                                    3.f * i);
  }

  SpacePointContainerConfig spConfig;
  spConfig.useDetailedDoubleMeasurementInfo = true;
  SpacePointContainerOptions spOptions;

  Adapter adapter(externalCollection);
  SpacePointContainer<Adapter, detail::RefHolder> spContainer(
      spConfig, spOptions, adapter);

  BOOST_CHECK_EQUAL(spContainer.size(), nExternalPoints);
  BOOST_CHECK_EQUAL(spContainer.size(), externalCollection.size());
  BOOST_CHECK_EQUAL(spContainer.end() - spContainer.begin(), nExternalPoints);
  BOOST_CHECK_EQUAL(std::distance(spContainer.begin(), spContainer.end()),
                    nExternalPoints);
  BOOST_CHECK_EQUAL(std::distance(std::ranges::begin(spContainer),
                                  std::ranges::end(spContainer)),
                    nExternalPoints);

  using proxy_t =
      SpacePointProxy<SpacePointContainer<Adapter, detail::RefHolder>>;
  static_assert(
      std::same_as<typename decltype(spContainer)::ValueType, SpacePoint>);
  static_assert(
      std::same_as<typename decltype(spContainer)::ProxyType, proxy_t>);
  static_assert(
      std::same_as<typename decltype(spContainer)::value_type, proxy_t>);
  static_assert(
      std::same_as<typename proxy_t::ContainerType, decltype(spContainer)>);
  static_assert(std::same_as<typename proxy_t::ValueType, SpacePoint>);

  using iterator_t = detail::ContainerIterator<decltype(spContainer), proxy_t&,
                                               std::size_t, false>;
  using const_iterator_t =
      detail::ContainerIterator<decltype(spContainer), const proxy_t&,
                                std::size_t, true>;
  static_assert(
      std::same_as<iterator_t, typename decltype(spContainer)::iterator>);
  static_assert(std::same_as<const_iterator_t,
                             typename decltype(spContainer)::const_iterator>);
  static_assert(
      std::same_as<iterator_t, decltype(std::ranges::begin(spContainer))>);
  static_assert(
      std::same_as<const_iterator_t, decltype(std::cbegin(spContainer))>);

  std::size_t n = 0ul;
  for (const proxy_t& proxy : spContainer) {
    float refX = 1.f * n;
    float refY = 1.5f * n;
    float refZ = 2.f * n;
    float refCovR = 2.5f * n;
    float refCovZ = 3.f * n;
    float refRadius = std::hypot(refX, refY);
    float refPhi = std::atan2(refY, refX);

    BOOST_CHECK_EQUAL(proxy.index(), n);
    BOOST_CHECK_EQUAL(proxy.x(), refX);
    BOOST_CHECK_EQUAL(proxy.y(), refY);
    BOOST_CHECK_EQUAL(proxy.z(), refZ);
    BOOST_CHECK_EQUAL(proxy.radius(), refRadius);
    BOOST_CHECK_EQUAL(proxy.phi(), refPhi);
    BOOST_CHECK_EQUAL(proxy.varianceR(), refCovR);
    BOOST_CHECK_EQUAL(proxy.varianceZ(), refCovZ);

    const Vector3& topStripVector = proxy.topStripVector();
    const Vector3& bottomStripVector = proxy.bottomStripVector();
    const Vector3& stripCenterDistance = proxy.stripCenterDistance();
    const Vector3& topStripCenterPosition = proxy.topStripCenterPosition();

    for (std::size_t i = 0; i < 3; ++i) {
      BOOST_CHECK_EQUAL(topStripVector[i], 0.);
      BOOST_CHECK_EQUAL(bottomStripVector[i], 0.);
      BOOST_CHECK_EQUAL(stripCenterDistance[i], 0.);
      BOOST_CHECK_EQUAL(topStripCenterPosition[i], 0.);
    }

    const SpacePoint& sp = proxy.externalSpacePoint();
    BOOST_CHECK_EQUAL(&sp, &externalCollection[n]);

    ++n;
  }
  BOOST_CHECK_EQUAL(n, nExternalPoints);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
