// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <chrono>
#include <exception>
#include <future>
#include <memory>
#include <stdexcept>
#include <thread>
#include <utility>

#include "NativeWorkerExchange.hpp"

using namespace ActsTests::NativeWorkers;
using namespace std::chrono_literals;

BOOST_AUTO_TEST_CASE(CancelResolvesQueuedAndRunningRequests) {
  auto exchange = std::make_shared<Exchange>(2);
  auto running = exchange->submit({11, 1, 0.});
  auto job = exchange->take();
  auto queued = exchange->submit({12, 2, 0.});
  BOOST_CHECK(running.wait_for(1ms) == std::future_status::timeout);
  BOOST_CHECK(queued.wait_for(1ms) == std::future_status::timeout);
  BOOST_CHECK_EQUAL(exchange->maxOutstanding(), 2u);

  std::promise<void> entering;
  auto started = entering.get_future();
  std::packaged_task<void()> submit(
      [exchange, entering = std::move(entering)]() mutable {
        entering.set_value();
        exchange->submit({13, 3, 0.});
      });
  auto admitted = submit.get_future();
  // The caller owns its exchange. A missing wakeup must fail by timeout
  // without blocking the test in a thread or future destructor.
  std::thread([submit = std::move(submit)]() mutable {
    submit.make_ready_at_thread_exit();
  }).detach();
  BOOST_REQUIRE(started.wait_for(5s) == std::future_status::ready);
  BOOST_CHECK(admitted.wait_for(50ms) == std::future_status::timeout);

  exchange->cancel(std::make_exception_ptr(std::domain_error("run aborted")));
  BOOST_REQUIRE(admitted.wait_for(5s) == std::future_status::ready);
  BOOST_CHECK_THROW(admitted.get(), std::domain_error);
  // The worker was already transporting this job when cancellation won.
  exchange->complete(job, {11, 0., {}});
  BOOST_REQUIRE(running.wait_for(0ms) == std::future_status::ready);
  BOOST_REQUIRE(queued.wait_for(0ms) == std::future_status::ready);
  BOOST_CHECK_THROW(running.get(), std::domain_error);
  BOOST_CHECK_THROW(queued.get(), std::domain_error);
  BOOST_CHECK(!exchange->take());
}

BOOST_AUTO_TEST_CASE(CloseDrainsQueuedAndRunningRequests) {
  Exchange exchange(2);
  auto first = exchange.submit({11, 1, 0.});
  auto running = exchange.take();
  auto second = exchange.submit({12, 2, 0.});
  exchange.close();
  BOOST_CHECK_THROW(exchange.submit({13, 3, 0.}), std::runtime_error);

  auto queued = exchange.take();
  BOOST_REQUIRE(queued);
  BOOST_CHECK_EQUAL(queued->request.eventId, 12);
  exchange.complete(queued, {12, 0., {}});
  BOOST_REQUIRE(second.wait_for(0ms) == std::future_status::ready);
  BOOST_CHECK_EQUAL(second.get().eventId, 12);
  BOOST_CHECK(!exchange.take());
  BOOST_CHECK(first.wait_for(1ms) == std::future_status::timeout);
  exchange.complete(running, {11, 0., {}});
  BOOST_REQUIRE(first.wait_for(0ms) == std::future_status::ready);
  BOOST_CHECK_EQUAL(first.get().eventId, 11);
}
