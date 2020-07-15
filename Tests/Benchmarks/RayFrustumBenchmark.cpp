// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
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
#include "Acts/Tests/CommonHelpers/BenchmarkTools.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts;

struct O {};
using Box = Acts::AxisAlignedBoundingBox<O, float, 3>;
using VertexType = Box::VertexType;
using vertex_array_type = Box::vertex_array_type;
using value_type = Box::value_type;
using Frustum3 = Frustum<float, 3, 4>;

bool optIntersect1(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();
  const vertex_array_type& idir = ray.idir();

  vertex_array_type t0s = (box.max() - origin).array() * idir;
  vertex_array_type t1s = (box.min() - origin).array() * idir;

  vertex_array_type tsmaller = t0s.min(t1s);
  vertex_array_type tbigger = t0s.max(t1s);

  value_type tmin = tsmaller.maxCoeff();
  value_type tmax = tbigger.minCoeff();

  return tmin < tmax && tmax > 0.0;
}

bool naiveIntersect1(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();

  const vertex_array_type& d = ray.dir();

  double tmin = -INFINITY, tmax = INFINITY;
  if (d.x() != 0.0) {
    double tx1 = (box.min().x() - origin.x()) / d.x();
    double tx2 = (box.max().x() - origin.x()) / d.x();

    tmin = std::max(tmin, std::min(tx1, tx2));
    tmax = std::min(tmax, std::max(tx1, tx2));
  }

  if (d.y() != 0.0) {
    double ty1 = (box.min().y() - origin.y()) / d.y();
    double ty2 = (box.max().y() - origin.y()) / d.y();

    tmin = std::max(tmin, std::min(ty1, ty2));
    tmax = std::min(tmax, std::max(ty1, ty2));
  }

  if (d.z() != 0.0) {
    double tz1 = (box.min().z() - origin.z()) / d.z();
    double tz2 = (box.max().z() - origin.z()) / d.z();

    tmin = std::max(tmin, std::min(tz1, tz2));
    tmax = std::min(tmax, std::max(tz1, tz2));
  }

  return tmax > tmin && tmax > 0.0;
}

bool naiveIntersect2(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();

  const vertex_array_type& d = ray.dir();

  double tmin = -INFINITY, tmax = INFINITY;

  for (size_t i = 0; i < 3; i++) {
    if (d[i] == 0.0)
      continue;
    double t1 = (box.min()[i] - origin[i]) / d[i];
    double t2 = (box.max()[i] - origin[i]) / d[i];
    tmin = std::max(tmin, std::min(t1, t2));
    tmax = std::min(tmax, std::max(t1, t2));
  }

  return tmax > tmin && tmax > 0.0;
}

bool naiveIntersect3(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();
  const vertex_array_type& d = ray.dir();

  double tx1 = (box.min().x() - origin.x()) / d.x();
  double tx2 = (box.max().x() - origin.x()) / d.x();
  double tmin = std::min(tx1, tx2);
  double tmax = std::max(tx1, tx2);

  double ty1 = (box.min().y() - origin.y()) / d.y();
  double ty2 = (box.max().y() - origin.y()) / d.y();
  tmin = std::max(tmin, std::min(ty1, ty2));
  tmax = std::min(tmax, std::max(ty1, ty2));

  double tz1 = (box.min().z() - origin.z()) / d.z();
  double tz2 = (box.max().z() - origin.z()) / d.z();
  tmin = std::max(tmin, std::min(tz1, tz2));
  tmax = std::min(tmax, std::max(tz1, tz2));

  return tmax > tmin && tmax > 0.0;
}

bool naiveIntersect4(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();
  const vertex_array_type& id = ray.idir();

  double tx1 = (box.min().x() - origin.x()) * id.x();
  double tx2 = (box.max().x() - origin.x()) * id.x();
  double tmin = std::min(tx1, tx2);
  double tmax = std::max(tx1, tx2);

  double ty1 = (box.min().y() - origin.y()) * id.y();
  double ty2 = (box.max().y() - origin.y()) * id.y();
  tmin = std::max(tmin, std::min(ty1, ty2));
  tmax = std::min(tmax, std::max(ty1, ty2));

  double tz1 = (box.min().z() - origin.z()) * id.z();
  double tz2 = (box.max().z() - origin.z()) * id.z();
  tmin = std::max(tmin, std::min(tz1, tz2));
  tmax = std::min(tmax, std::max(tz1, tz2));

  return tmax > tmin && tmax > 0.0;
}

bool naiveIntersect5(const Box& box, const Ray<float, 3>& ray) {
  const VertexType& origin = ray.origin();
  const vertex_array_type& id = ray.idir();
  double tmin = -INFINITY, tmax = INFINITY;

  for (size_t i = 0; i < 3; i++) {
    double t1 = (box.min()[i] - origin[i]) * id[i];
    double t2 = (box.max()[i] - origin[i]) * id[i];
    tmin = std::max(tmin, std::min(t1, t2));
    tmax = std::min(tmax, std::max(t1, t2));
  }

  return tmax > tmin && tmax > 0.0;
}

bool frustOpt1(const Box& box, const Frustum3& fr) {
  const auto& normals = fr.normals();
  const vertex_array_type fr_vmin = box.min() - fr.origin();
  const vertex_array_type fr_vmax = box.max() - fr.origin();

  VertexType p_vtx;
  for (size_t i = 0; i < fr.sides + 1; i++) {
    const VertexType& normal = normals[i];

    p_vtx = (normal.array() < 0).template cast<value_type>() * fr_vmin +
            (normal.array() >= 0).template cast<value_type>() * fr_vmax;

    if (p_vtx.dot(normal) < 0) {
      return false;
    }
  }
  return true;
}

template <typename frustum_t>
bool frustOpt2(const Box& box, const frustum_t& fr) {
  const auto& normals = fr.normals();
  const vertex_array_type fr_vmin = box.min() - fr.origin();
  const vertex_array_type fr_vmax = box.max() - fr.origin();

  auto calc = [&](const auto& normal) {
    return (normal.array() < 0).template cast<value_type>() * fr_vmin +
           (normal.array() >= 0).template cast<value_type>() * fr_vmax;
  };

  VertexType p_vtx;

  p_vtx = calc(normals[0]);
  if (p_vtx.dot(normals[0]) < 0) {
    return false;
  }

  p_vtx = calc(normals[1]);
  if (p_vtx.dot(normals[1]) < 0) {
    return false;
  }

  p_vtx = calc(normals[2]);
  if (p_vtx.dot(normals[2]) < 0) {
    return false;
  }

  if constexpr (fr.sides > 2) {
    p_vtx = calc(normals[3]);
    if (p_vtx.dot(normals[3]) < 0) {
      return false;
    }
  }

  if constexpr (fr.sides > 3) {
    p_vtx = calc(normals[4]);
    if (p_vtx.dot(normals[4]) < 0) {
      return false;
    }
  }

  if constexpr (fr.sides > 4) {
    for (size_t i = 5; i <= fr.sides; i++) {
      const VertexType& normal = normals[i];

      p_vtx = calc(normal);
      if (p_vtx.dot(normal) < 0) {
        return false;
      }
    }
  }

  return true;
}

bool frustOpt3(const Box& box, const Frustum3& fr) {
  const auto& normals = fr.normals();
  const vertex_array_type fr_vmin = box.min() - fr.origin();
  const vertex_array_type fr_vmax = box.max() - fr.origin();

  VertexType p_vtx;
  bool result = true;
  for (size_t i = 0; i < fr.sides + 1; i++) {
    const VertexType& normal = normals[i];

    p_vtx = (normal.array() < 0).template cast<value_type>() * fr_vmin +
            (normal.array() >= 0).template cast<value_type>() * fr_vmax;

    result = result && (p_vtx.dot(normal) >= 0);
  }
  return result;
}

template <typename frustum_t>
bool frustOpt4(const Box& box, const frustum_t& fr) {
  const auto& normals = fr.normals();
  const vertex_array_type fr_vmin = box.min() - fr.origin();
  const vertex_array_type fr_vmax = box.max() - fr.origin();

  auto calc = [&](const auto& normal) {
    return (normal.array() < 0).template cast<value_type>() * fr_vmin +
           (normal.array() >= 0).template cast<value_type>() * fr_vmax;
  };

  VertexType p_vtx;
  bool result = true;

  if constexpr (fr.sides > 4) {
    for (size_t i = 0; i < fr.sides + 1; i++) {
      const VertexType& normal = normals[i];

      p_vtx = calc(normal);
      result = result && (p_vtx.dot(normal) >= 0);
    }
  } else {
    p_vtx = calc(normals[0]);
    result = result && (p_vtx.dot(normals[0]) >= 0);

    p_vtx = calc(normals[1]);
    result = result && (p_vtx.dot(normals[1]) >= 0);

    p_vtx = calc(normals[2]);
    result = result && (p_vtx.dot(normals[2]) >= 0);

    if constexpr (fr.sides > 2) {
      p_vtx = calc(normals[3]);
      result = result && (p_vtx.dot(normals[3]) >= 0);

      if constexpr (fr.sides > 3) {
        p_vtx = calc(normals[4]);
        result = result && (p_vtx.dot(normals[4]) >= 0);
      }
    }
  }

  return result;
}

int main(int /*argc*/, char** /*argv[]*/) {
  size_t n = 10000;

  std::vector<Frustum3> frustums{n, Frustum3{{0, 0, 0}, {1, 0, 0}, M_PI / 2.}};
  std::mt19937 rng(42);
  std::uniform_real_distribution<float> dir(0, 1);
  std::uniform_real_distribution<float> loc(-10, 10);
  std::uniform_real_distribution<float> ang(M_PI / 10., M_PI / 4.);

  std::generate(frustums.begin(), frustums.end(), [&]() {
    const Vector3F d{dir(rng), dir(rng), dir(rng)};
    const Vector3F l{loc(rng), loc(rng), loc(rng)};
    return Frustum3{l, d.normalized(), ang(rng)};
  });

  using Ray3 = Ray<float, 3>;
  std::vector<Ray3> rays{n, Ray3{Vector3F{0, 0, 0}, Vector3F{1, 0, 0}}};
  std::generate(rays.begin(), rays.end(), [&]() {
    const Vector3F d{dir(rng), dir(rng), dir(rng)};
    const Vector3F l{loc(rng), loc(rng), loc(rng)};
    return Ray3{l, d.normalized()};
  });

  Box box{nullptr, {0, 0, 0}, Box::Size{{1, 2, 3}}};

  std::cout << "Make sure ray implementations are identical" << std::endl;
  for (const auto& ray : rays) {
    std::vector<bool> results{
        optIntersect1(box, ray),   naiveIntersect1(box, ray),
        naiveIntersect2(box, ray), naiveIntersect3(box, ray),
        naiveIntersect4(box, ray), naiveIntersect5(box, ray),
    };
    bool all =
        std::all_of(results.begin(), results.end(), [](bool r) { return r; });
    bool none =
        std::none_of(results.begin(), results.end(), [](bool r) { return r; });

    if (!all && !none) {
      std::cerr << "Discrepancy: " << std::endl;
      for (size_t i = 0; i < results.size(); i++) {
        std::cerr << " - #" << i << ": " << results.at(i) << std::endl;
      }

      box.toStream(std::cerr);
      std::cerr << std::endl;
      std::cerr << "Ray: [" << ray.origin().transpose() << "], ["
                << ray.dir().transpose() << "]" << std::endl;
      return -1;
    }
  }
  std::cout << "Seems ok" << std::endl;

  std::cout << "Benchmarking rays nominal: " << std::flush;
  auto bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return box.intersect(ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Benchmarking rays naive 1: " << std::flush;
  bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return naiveIntersect1(box, ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Benchmarking rays naive 2: " << std::flush;
  bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return naiveIntersect2(box, ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Benchmarking rays naive 3: " << std::flush;
  bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return naiveIntersect3(box, ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Benchmarking rays naive 4: " << std::flush;
  bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return naiveIntersect4(box, ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Benchmarking rays naive 5: " << std::flush;
  bench_result = Acts::Test::microBenchmark(
      [&](const auto& ray) { return naiveIntersect5(box, ray); }, rays);
  std::cout << bench_result << std::endl;

  std::cout << "Make sure frustum implementations are identical" << std::endl;
  for (const auto& fr : frustums) {
    std::vector<bool> results{
        box.intersect(fr),  frustOpt1(box, fr), frustOpt2(box, fr),
        frustOpt3(box, fr), frustOpt4(box, fr),
    };
    bool all =
        std::all_of(results.begin(), results.end(), [](bool r) { return r; });
    bool none =
        std::none_of(results.begin(), results.end(), [](bool r) { return r; });

    if (!all && !none) {
      std::cerr << "Discrepancy: " << std::endl;
      for (size_t i = 0; i < results.size(); i++) {
        std::cerr << " - #" << i << ": " << results.at(i) << std::endl;
      }

      box.toStream(std::cerr);
      std::cerr << std::endl;
      std::cerr << "Ray: [" << fr.origin().transpose() << "], ["
                << fr.dir().transpose() << "]" << std::endl;
      return -1;
    }
  }
  std::cout << "Seems ok" << std::endl;

  std::cout << "Benchmarking frust nominal: " << std::flush;
  auto fr_bench_result = Acts::Test::microBenchmark(
      [&](const auto& fr) { return box.intersect(fr); }, frustums);
  std::cout << fr_bench_result << std::endl;

  std::cout << "Benchmarking frust opt 1: " << std::flush;
  fr_bench_result = Acts::Test::microBenchmark(
      [&](const auto& fr) { return frustOpt1(box, fr); }, frustums);
  std::cout << fr_bench_result << std::endl;

  std::cout << "Benchmarking frust opt 2: " << std::flush;
  fr_bench_result = Acts::Test::microBenchmark(
      [&](const auto& fr) { return frustOpt2(box, fr); }, frustums);
  std::cout << fr_bench_result << std::endl;

  std::cout << "Benchmarking frust opt 3: " << std::flush;
  fr_bench_result = Acts::Test::microBenchmark(
      [&](const auto& fr) { return frustOpt3(box, fr); }, frustums);
  std::cout << fr_bench_result << std::endl;

  std::cout << "Benchmarking frust opt 4: " << std::flush;
  fr_bench_result = Acts::Test::microBenchmark(
      [&](const auto& fr) { return frustOpt4(box, fr); }, frustums);
  std::cout << fr_bench_result << std::endl;

  return 0;
}
