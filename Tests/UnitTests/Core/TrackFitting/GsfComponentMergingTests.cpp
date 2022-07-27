// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <random>

#include <Eigen/Eigenvalues> 

#include <boost/test/unit_test.hpp>

#include "Acts/TrackFitting/detail/KLMixtureReduction.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"

#define CHECK_CLOSE_MATRIX(a, b, t) BOOST_CHECK(((a-b).array().abs() < t).all())

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::UnitConstants;

struct DummyComponent {
  Acts::ActsScalar weight;
  Acts::BoundVector boundPars;
  std::optional<Acts::BoundSymMatrix> boundCov;
};

template<typename T, int D>
class MultivariateNormalDistribution
{
public:
    using Vector = Eigen::Matrix<T, D, 1>;
    using Matrix = Eigen::Matrix<T, D, D>;
  
private:
    Vector m_mean;
    Matrix m_transform;
    
public:
    MultivariateNormalDistribution(Vector const& mean, Matrix const& cov)
        : m_mean(mean)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(cov);
        m_transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    template<typename generator_t>
    Vector operator()(generator_t &gen) const
    {
      std::normal_distribution<T> normal;
        return m_mean + m_transform * Vector{}.unaryExpr([&](auto) { return normal(gen); });
    }
};

auto meanCovFromData(const std::vector<DummyComponent> &cmps, std::size_t n_samples, std::mt19937 &gen) {
  using MultiNormal = MultivariateNormalDistribution<double, eBoundSize>;
  
  std::vector<std::pair<double, MultiNormal>> dists;
  for(const auto &cmp : cmps) {
    dists.emplace_back(cmp.weight, MultiNormal(cmp.boundPars, *cmp.boundCov));
  }
  
  auto sample = [&](){
    MultiNormal::Vector ret = MultiNormal::Vector::Zero();
    
    for(const auto &[w, dist] : dists) {
      ret += w * dist(gen);
    }
    
    return ret;
  };
  
  std::vector<BoundVector> samples(n_samples);
  std::generate(samples.begin(), samples.end(), sample);
  
  BoundVector mean = BoundVector::Zero();
  for(const auto &x : samples) {
    mean += x;
  }
  mean /= samples.size();
  
  BoundSymMatrix cov = BoundSymMatrix::Zero();
  for(const auto &x : samples) {
    cov += (x - mean) * (x - mean).transpose();
  }
  cov /= samples.size();
  
  return std::make_tuple(mean, cov);
}

BoundVector meanFromFree(const std::vector<DummyComponent> &cmps, const Surface &surface) {
  
  if( surface.type() == Surface::Cylinder ) {
    throw std::runtime_error("Cylinder surface not yet supported");
  }
  
  if( surface.type() == Surface::Cone ) {
    throw std::runtime_error("Cone surface not supported");
  }
  
  FreeVector mean = FreeVector::Zero();
  
  for(const auto &cmp : cmps) {
    mean += cmp.weight * detail::transformBoundToFreeParameters(surface, GeometryContext{}, cmp.boundPars);
  }
  
  mean.segment<3>(eFreeDir0).normalize();
  
  return *detail::transformFreeToBoundParameters(mean, surface, GeometryContext{});
}



BOOST_AUTO_TEST_CASE(test_plane_surface) {
  const auto desc = detail::AngleDescription::Default{};
  const auto proj = Identity{};
  
  Acts::BoundSymMatrix cov_ref = Acts::BoundSymMatrix::Zero();
  cov_ref.diagonal() = BoundVector::Ones() / 10;
  cov_ref.diagonal(+1) = ActsVector<5>::Ones() / 100;
  cov_ref.diagonal(-1) = ActsVector<5>::Ones() / 100;
  
  const auto surface = Surface::makeShared<PlaneSurface>(Vector3{0,0,0}, Vector3{1,0,0});
  
  for(auto phi : { -180_degree, 0_degree, 180_degree }) {
    for(auto theta : { 5_degree, 90_degree, 175_degree }) {
      // Go create mixture with 4 cmps
      std::vector<DummyComponent> cmps;
      
      for(auto dphi : { -10_degree, 10_degree }) {
        for(auto dtheta : {-10_degree, 10_degree }) {
          DummyComponent a;
          a.weight = 1./4.;
          a.boundPars = BoundVector::Random(); // between -1 and 1
          a.boundPars[eBoundPhi] = phi + dphi;
          a.boundPars[eBoundTheta] = theta + dtheta;         
          a.boundCov = cov_ref;
          
          cmps.push_back(a);
        }
      }
      
      const auto [mean_test, cov_test] = detail::combineBoundGaussianMixture(cmps.begin(), cmps.end(), proj, desc);
      
//       const auto mean_ref = meanFromFree(cmps, *surface);
      
      std::mt19937 gen(std::random_device{}());
      
      const auto [mean_ref, cov_data] = meanCovFromData(cmps, 1'000'000, gen);
      
//       std::cout << "phi: " << phi / degree << ", theta: " << theta / degree << "\n";
//       std::cout << (mean_test.transpose() / degree) << "\n";
//       std::cout << (mean_ref.transpose() / degree) << "\n";
//       std::cout << (mean_ref - mean_test).transpose() / degree << "\n";
//       std::cout << (mean_ref.array() / mean_test.array()).transpose() << "\n";
      
      std::cout << "cov_test:\n" << *cov_test << "\n";
      std::cout << "cov_ref:\n" << cov_data << "\n";
      std::cout << "diff:\n" << (*cov_test - cov_data).array().abs() << "\n";
      
      CHECK_CLOSE_MATRIX(mean_test, mean_ref, 1.e-1_degree);
      CHECK_CLOSE_MATRIX(*cov_test, cov_data, 1.e-3);
      return;
    }
  }
}

/*
BOOST_AUTO_TEST_CASE(test_merge_two_equal_components) {
  DummyComponent a;
  a.boundPars = Acts::BoundVector::Random();
  a.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
  *a.boundCov *= a.boundCov->transpose();
  a.weight = 0.5;

  DummyComponent c = Acts::detail::mergeComponents(a, a, Identity{});
  BOOST_CHECK(c.boundPars == a.boundPars);
  BOOST_CHECK(*c.boundCov == *a.boundCov);
  BOOST_CHECK(c.weight == 1.0);
}

BOOST_AUTO_TEST_CASE(test_merge_two_different_components) {
  DummyComponent a;
  a.boundPars = Acts::BoundVector::Random();
  a.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
  *a.boundCov *= a.boundCov->transpose();
  a.weight = 0.5;

  DummyComponent b;
  b.boundPars = Acts::BoundVector::Random();
  b.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
  *b.boundCov *= b.boundCov->transpose();
  b.weight = 0.5;

  DummyComponent c = Acts::detail::mergeComponents(a, b, Identity{});
  BOOST_CHECK(c.boundPars == 0.5 * (a.boundPars + b.boundPars));
  BOOST_CHECK(c.weight == 1.0);
}

BOOST_AUTO_TEST_CASE(test_component_reduction_equal) {
  const std::size_t NCompsBefore = 10;

  // Create start state
  std::vector<DummyComponent> cmps;

  DummyComponent a;
  a.boundPars = Acts::BoundVector::Random();
  a.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
  *a.boundCov *= a.boundCov->transpose();
  a.weight = 1.0 / NCompsBefore;

  for (auto i = 0ul; i < NCompsBefore; ++i) {
    cmps.push_back(a);
  }

  const double weightSumBefore = std::accumulate(
      cmps.begin(), cmps.end(), 0.0,
      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  BOOST_CHECK_CLOSE(weightSumBefore, 1.0, 0.0001);

  // Combine
  while (cmps.size() >= 2) {
    auto merge_iter_a = cmps.begin();
    auto merge_iter_b = std::next(cmps.begin());

    *merge_iter_a =
        Acts::detail::mergeComponents(*merge_iter_a, *merge_iter_b, Identity{});
    cmps.erase(merge_iter_b);

    const auto mean = std::accumulate(
        cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
        [](auto sum, const auto &cmp) -> Acts::BoundVector {
          return sum + cmp.weight * cmp.boundPars;
        });

    const double weightSum = std::accumulate(
        cmps.begin(), cmps.end(), 0.0,
        [](auto sum, const auto &cmp) { return sum + cmp.weight; });

    BOOST_CHECK((mean - a.boundPars).cwiseAbs().all() < 1.e-4);
    BOOST_CHECK_CLOSE(weightSum, 1.0, 0.0001);

    if (cmps.size() == 1) {
      BOOST_CHECK_CLOSE(weightSum, merge_iter_a->weight, 0.0001);
      BOOST_CHECK((a.boundPars - merge_iter_a->boundPars).cwiseAbs().all() <
                  1.e-4);
      BOOST_CHECK((*a.boundCov - *merge_iter_a->boundCov).cwiseAbs().all() <
                  1.e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(test_component_reduction_different) {
  const std::size_t NCompsBefore = 10;

  // Create start state
  std::vector<DummyComponent> cmps;

  for (auto i = 0ul; i < NCompsBefore; ++i) {
    DummyComponent a;
    a.boundPars = Acts::BoundVector::Random();
    a.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
    *a.boundCov *= a.boundCov->transpose();
    a.weight = 1.0 / NCompsBefore;
    cmps.push_back(a);
  }

  // Determine mean
  const auto meanBefore = std::accumulate(
      cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
      [](auto sum, const auto &cmp) -> Acts::BoundVector {
        return sum + cmp.weight * cmp.boundPars;
      });

  const double weightSumBefore = std::accumulate(
      cmps.begin(), cmps.end(), 0.0,
      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  BOOST_CHECK_CLOSE(weightSumBefore, 1.0, 0.0001);

  // Combine
  while (cmps.size() >= 2) {
    auto merge_iter_a = cmps.begin();
    auto merge_iter_b = std::next(cmps.begin());

    *merge_iter_a =
        Acts::detail::mergeComponents(*merge_iter_a, *merge_iter_b, Identity{});
    cmps.erase(merge_iter_b);

    const auto mean = std::accumulate(
        cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
        [](auto sum, const auto &cmp) -> Acts::BoundVector {
          return sum + cmp.weight * cmp.boundPars;
        });

    const double weightSum = std::accumulate(
        cmps.begin(), cmps.end(), 0.0,
        [](auto sum, const auto &cmp) { return sum + cmp.weight; });

    BOOST_CHECK((mean - meanBefore).cwiseAbs().all() < 1.e-4);
    BOOST_CHECK_CLOSE(weightSum, 1.0, 0.0001);
  }
}

BOOST_AUTO_TEST_CASE(test_kl_mixture_reduction) {
  const std::size_t NCompsBefore = 10;
  const std::size_t NCompsAfter = 5;

  // Create start state
  std::vector<DummyComponent> cmps;

  for (auto i = 0ul; i < NCompsBefore; ++i) {
    DummyComponent a;
    a.boundPars = Acts::BoundVector::Random();
    a.boundCov = Acts::BoundSymMatrix::Random().cwiseAbs();
    *a.boundCov *= a.boundCov->transpose();
    a.weight = 1.0 / NCompsBefore;
    cmps.push_back(a);
  }

  // Determine mean
  const auto meanBefore = std::accumulate(
      cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
      [](auto sum, const auto &cmp) -> Acts::BoundVector {
        return sum + cmp.weight * cmp.boundPars;
      });

  const double weightSumBefore = std::accumulate(
      cmps.begin(), cmps.end(), 0.0,
      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  BOOST_CHECK_CLOSE(weightSumBefore, 1.0, 0.0001);

  Acts::detail::reduceWithKLDistance(cmps, NCompsAfter, Identity{});

  const auto meanAfter = std::accumulate(
      cmps.begin(), cmps.end(), Acts::BoundVector::Zero().eval(),
      [](auto sum, const auto &cmp) -> Acts::BoundVector {
        return sum + cmp.weight * cmp.boundPars;
      });

  const double weightSumAfter = std::accumulate(
      cmps.begin(), cmps.end(), 0.0,
      [](auto sum, const auto &cmp) { return sum + cmp.weight; });

  BOOST_CHECK_CLOSE(weightSumAfter, 1.0, 0.0001);
  BOOST_CHECK((meanAfter - meanBefore).cwiseAbs().all() < 1.e-4);
  BOOST_CHECK(cmps.size() == NCompsAfter);
}
*/
