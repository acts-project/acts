// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts;

constexpr int NTESTS = 2 << 15;
// trapezoidal area
static const Vector2D POLY[]
    = {{0.4, 0.25}, {0.6, 0.25}, {0.8, 0.75}, {0.2, 0.75}};

std::vector<Vector2D>
make_random_points()
{
  std::mt19937                           rng(42);
  std::uniform_real_distribution<double> axis(0, 1);

  auto rnd_axis  = std::bind(axis, rng);
  auto rnd_point = [&]() { return Vector2D(rnd_axis(), rnd_axis()); };

  std::vector<Vector2D> points(NTESTS);
  std::generate(points.begin(), points.end(), rnd_point);
  return points;
}

// inline bool
// isInsideOld(const Vector2D& point, const BoundaryCheck& check)
//{
//  // reproduce most of the steps e.g. in TrapezoidBounds::inside

//  // rotate to a vector normal to the input
//  ActsMatrixD<2, 2> toNormal;
//  toNormal << 0, -1, 1, 0;

//  // compute KDOP and axes for surface polygon
//  std::vector<Vector2D> points = {POLY[0], POLY[1], POLY[2], POLY[3]};
//  std::vector<Vector2D> axes   = {toNormal * (POLY[1] - POLY[0]),
//                                toNormal * (POLY[2] - POLY[1]),
//                                toNormal * (POLY[3] - POLY[2])};
//  std::vector<KDOP> limits(3);
//  std::vector<KDOP> limitsErr(3);

//  check.ComputeKDOP(points, axes, limits);
//  check.ComputeKDOP(check.EllipseToPoly(3), axes, limitsErr);
//  // check if KDOPs overlap and return result
//  return check.TestKDOPKDOP(limits, limitsErr);
//}

inline bool
isInside(const Vector2D& point, const BoundaryCheck& check)
{
  return check.isInside(point, POLY);
}

struct StopWatch
{
  using Clock     = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;

  TimePoint start;

  StopWatch() : start(Clock::now()) {}

  void
  finish(size_t n_trials, const char* name)
  {
    auto stop = Clock::now();
    std::chrono::duration<double, std::micro> total = stop - start;
    std::chrono::duration<double, std::micro> perTrial
        = (stop - start) / n_trials;

    std::cout << name << ":\n";
    std::cout << "  trials: " << n_trials << '\n';
    std::cout << "  time total: " << total.count() << " us\n";
    std::cout << "  time per test: " << perTrial.count() << " us\n";
  }
};

int
main(int /*argc*/, char** /*argv[]*/)
{
  using std::cout;

  // absolute check w/o tolerance
  {
    auto          points = make_random_points();
    BoundaryCheck check(true);
    int           n_inside = 0;
    StopWatch     watch;
    for (const auto& point : points) {
      n_inside += (isInside(point, check) ? 1 : 0);
    }
    watch.finish(points.size(), "absolute check w/o tolerance");
  }
  // check w/ covariance
  {
    auto              points = make_random_points();
    ActsSymMatrixD<2> cov;
    cov << 0.2, 0.02, 0.15, 0.02;
    BoundaryCheck check(cov, 3.0);  // 3-sigma cut
    int           n_inside = 0;
    StopWatch     watch;
    for (const auto& point : points) {
      n_inside += (isInside(point, check) ? 1 : 0);
    }
    watch.finish(points.size(), "check w/ covariance");
  }
  // check w/ covariance
  //  {
  //    auto              points = make_random_points();
  //    ActsSymMatrixD<2> cov;
  //    cov << 0.2, 0.02, 0.15, 0.02;
  //    BoundaryCheck check(cov, 3.0);  // 3-sigma cut
  //    int           n_inside = 0;
  //    StopWatch     watch;
  //    for (const auto& point : points) {
  //      n_inside += (isInsideOld(point, check) ? 1 : 0);
  //    }
  //    watch.finish(points.size(), "old check w/ covariance");
  //  }

  return 0;
}
