// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/detail/WorkerPool.hpp"

#include <algorithm>
#include <atomic>
#include <iterator>
#include <stdexcept>
#include <thread>
#include <vector>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(WorkerPoolVisitsEveryIndexOnce) {
  for (std::size_t nThreads : {0u, 1u, 2u, 4u, 7u}) {
    detail::WorkerPool pool(nThreads);
    BOOST_CHECK_EQUAL(pool.size(), std::max<std::size_t>(nThreads, 1));
    for (std::size_t n : {0u, 1u, 2u, 13u, 1000u}) {
      for (std::size_t grainSize : {0u, 1u, 3u, 5000u}) {
        std::vector<std::atomic<int>> visits(n);
        pool.parallelFor(
            n, [&](std::size_t i) { visits[i].fetch_add(1); }, grainSize);
        for (const auto& v : visits) {
          BOOST_CHECK_EQUAL(v.load(), 1);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(WorkerPoolUsesThreads) {
  detail::WorkerPool pool(4);
  std::vector<std::thread::id> ids(64);
  std::atomic<int> waiting{0};
  pool.parallelFor(
      ids.size(),
      [&](std::size_t i) {
        ids[i] = std::this_thread::get_id();
        // Block the first indices until every thread has claimed one, so
        // the loop cannot be drained by a single thread
        if (i < 4) {
          waiting.fetch_add(1);
          while (waiting.load() < 4) {
            std::this_thread::yield();
          }
        }
      },
      1);
  std::ranges::sort(ids);
  auto unique = std::ranges::unique(ids);
  BOOST_CHECK_EQUAL(std::distance(ids.begin(), unique.begin()), 4);
}

BOOST_AUTO_TEST_CASE(WorkerPoolRethrows) {
  detail::WorkerPool pool(4);
  std::atomic<int> calls{0};
  BOOST_CHECK_THROW(pool.parallelFor(1000,
                                     [&](std::size_t i) {
                                       calls.fetch_add(1);
                                       if (i == 17) {
                                         throw std::runtime_error("failure");
                                       }
                                     }),
                    std::runtime_error);
  BOOST_CHECK_GE(calls.load(), 1);

  // The pool is still usable afterwards
  std::atomic<std::size_t> sum{0};
  pool.parallelFor(100, [&](std::size_t i) { sum.fetch_add(i); });
  BOOST_CHECK_EQUAL(sum.load(), 4950u);
}

BOOST_AUTO_TEST_CASE(WorkerPoolNestedAndConcurrentCalls) {
  detail::WorkerPool pool(4);

  // A nested call from a loop body runs serially instead of dead-locking
  std::atomic<std::size_t> nested{0};
  pool.parallelFor(8, [&](std::size_t) {
    pool.parallelFor(10, [&](std::size_t) { nested.fetch_add(1); });
  });
  BOOST_CHECK_EQUAL(nested.load(), 80u);

  // Cross-pool nesting must not lose the outer pool's nesting guard, or
  // attempt to reacquire its call mutex from either pool's loop body.
  detail::WorkerPool other(2);
  nested.store(0);
  pool.parallelFor(8, [&](std::size_t) {
    other.parallelFor(2, [&](std::size_t) {
      pool.parallelFor(3, [&](std::size_t) { nested.fetch_add(1); });
    });
    pool.parallelFor(4, [&](std::size_t) { nested.fetch_add(1); });
  });
  BOOST_CHECK_EQUAL(nested.load(), 80u);

  // Calls from different threads are serialized
  std::atomic<std::size_t> total{0};
  std::vector<std::thread> callers;
  for (int c = 0; c < 3; ++c) {
    callers.emplace_back([&] {
      for (int r = 0; r < 20; ++r) {
        pool.parallelFor(100, [&](std::size_t) { total.fetch_add(1); });
      }
    });
  }
  for (auto& caller : callers) {
    caller.join();
  }
  BOOST_CHECK_EQUAL(total.load(), 3u * 20u * 100u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
