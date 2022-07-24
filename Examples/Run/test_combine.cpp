#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/detail/gaussian_mixture_helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#define PHI_THETA(dir) "(" << VectorHelpers::phi(dir) / u::degree << ", " << VectorHelpers::theta(dir) / u::degree << ")"

using namespace Acts::detail;
using namespace Acts;

namespace u = Acts::UnitConstants;

template <typename component_iterator_t, typename projector_t = Identity,
          typename angle_desc_t = AngleDescription::Default>
auto myCombineBoundGaussianMixture(
    const component_iterator_t begin, const component_iterator_t end,
    projector_t &&projector = projector_t{},
    const angle_desc_t &angle_desc = angle_desc_t{}) {
  throw_assert(std::distance(begin, end) > 0, "empty component range");

  using ret_type = std::tuple<BoundVector, std::optional<BoundSymMatrix>>;
  BoundVector mean = BoundVector::Zero();
  BoundSymMatrix cov1 = BoundSymMatrix::Zero();
  BoundSymMatrix cov2 = BoundSymMatrix::Zero();
  double sumOfWeights{0.0};

  const auto &[begin_weight, begin_pars, begin_cov] = projector(*begin);

  if (std::distance(begin, end) == 1) {
    return ret_type{begin_pars, *begin_cov};
  }
  
  Vector3 dir = Vector3::Zero();

  // clang-format off
  // x = \sum_{l} w_l * x_l
  // C = \sum_{l} w_l * C_l + \sum_{l} \sum_{m>l} w_l * w_m * (x_l - x_m)(x_l - x_m)^T
  // clang-format on
  for (auto l = begin; l != end; ++l) {
    const auto &[weight_l, pars_l, cov_l] = projector(*l);
    throw_assert(cov_l, "we require a covariance here");

    sumOfWeights += weight_l;
    mean += weight_l * pars_l;
    cov1 += weight_l * *cov_l;
    
    dir += weight_l * Acts::makeDirectionUnitFromPhiTheta(pars_l[eBoundPhi],
                                                          pars_l[eBoundTheta]);

    // Avoid problems with cyclic coordinates. The indices for these are taken
    // from the angle_description_t template parameter, and applied with the
    // following lambda.
    auto handleCyclicCoor = [&ref = begin_pars, &pars = pars_l,
                             &weight = weight_l, &mean = mean](auto desc) {
      const auto delta = (ref[desc.idx] - pars[desc.idx]) / desc.constant;
      
//       const auto degree = 360. / (2*M_PI);
//       
//       std::cout << "idx:  " << desc.idx << "\n";
//       std::cout << "ref:  " << ref[desc.idx] * degree << "\n";
//       std::cout << "pars:  " << pars[desc.idx] * degree << "\n";

      if (delta > M_PI) {
        mean[desc.idx] += (2 * M_PI) * weight * desc.constant;
      } else if (delta < -M_PI) {
        mean[desc.idx] -= (2 * M_PI) * weight * desc.constant;
      }
      
//       std::cout << "mean: " << mean[desc.idx] * degree << "\n---\n";
    };

    std::apply([&](auto... idx_desc) { (handleCyclicCoor(idx_desc), ...); },
               angle_desc);

    // For covariance we must loop over all other following components
    for (auto m = std::next(l); m != end; ++m) {
      const auto &[weight_m, pars_m, cov_m] = projector(*m);
      throw_assert(cov_m, "we require a covariance here");

      const BoundVector diff = pars_l - pars_m;
      cov2 += weight_l * weight_m * diff * diff.transpose();
    }
  }
  
  mean[eBoundTheta] = VectorHelpers::theta(dir);
  mean[eBoundPhi] = VectorHelpers::phi(dir);

  return ret_type{mean / sumOfWeights,
                  cov1 / sumOfWeights + cov2 / (sumOfWeights * sumOfWeights)};
}

int main() {
//     std::srand((unsigned int) std::time(0));
    std::srand(42);

    const auto scale = 1.;
    const auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3{0,0,0}, Acts::Vector3{1,0,0});

//     const auto scale = 100.;
//     const auto surface = Acts::Surface::makeShared<Acts::DiscSurface>(Acts::Transform3::Identity(), 0.0, scale);
    const Acts::GeometryContext gctx{};
    std::cout << std::tie(*surface, gctx) << "\n";
    int n=0;
    const int N=100;
    for(int i=0; i<N; ++i) {
        
        auto bound_pos = Acts::Vector2::Random().eval();
        bound_pos[0] *= scale;
        bound_pos[1] *= 2 * M_PI;

        Acts::FreeVector pars1 = Acts::FreeVector::Zero();
        Acts::FreeVector pars2 = Acts::FreeVector::Zero();

       pars1.segment<3>(Acts::eFreeDir0) = Acts::Vector3::Random().eval().normalized();
       pars2.segment<3>(Acts::eFreeDir0) = Acts::Vector3::Random().eval().normalized();
//          pars1.segment<3>(Acts::eFreeDir0) = Acts::Vector3{1,-1,1}.normalized();
//          pars2.segment<3>(Acts::eFreeDir0) = Acts::Vector3{1,1,1}.normalized();
                
//         const auto free_pos = surface->localToGlobal(gctx, bound_pos, pars1.segment<3>(Acts::eFreeDir0));
//         pars1.segment<3>(Acts::eFreePos0) = free_pos;
//         pars2.segment<3>(Acts::eFreePos0) = free_pos;

        const auto bound1 = *Acts::detail::transformFreeToBoundParameters(pars1, *surface, gctx);
        const auto bound2 = *Acts::detail::transformFreeToBoundParameters(pars2, *surface, gctx);
        
        std::array<std::tuple<double, Acts::BoundVector, std::optional<Acts::BoundSymMatrix>>, 2> a;

        a[0] = {0.5, bound1, Acts::BoundSymMatrix::Identity().eval() };
        a[1] = {0.5, bound2, Acts::BoundSymMatrix::Identity().eval() };

        const auto [res_bound, cov] = myCombineBoundGaussianMixture(a.begin(), a.end());

        const auto dir_from_bound = Acts::detail::transformBoundToFreeParameters(*surface, gctx, res_bound).segment<3>(Acts::eFreeDir0).normalized();
        const auto dir_from_free = (0.5 * pars1.segment<3>(Acts::eFreeDir0) + 0.5 * pars2.segment<3>(Acts::eFreeDir0)).normalized();

        const double res = (dir_from_bound.normalized() - dir_from_free.normalized()).array().abs().sum();

        if(res > 1.e-4 ) {
            std::cout << "bound1:      " << bound1.transpose() / u::degree << "\n";
            std::cout << "bound2:      " << bound2.transpose() / u::degree << "\n";
            std::cout << "bound_after: " << res_bound.transpose() / u::degree << "\n";
            std::cout << PHI_THETA(pars1.segment<3>(Acts::eFreeDir0)) << " + " 
                      << PHI_THETA(pars2.segment<3>(Acts::eFreeDir0)) << " -> " 
                      << PHI_THETA(dir_from_free) << "\n";
            std::cout << "res: " << res << "\n";
            std::cout << "dir1 before:    " << pars1.segment<3>(Acts::eFreeDir0).transpose() << "\n";
            std::cout << "dir2 before:    " << pars2.segment<3>(Acts::eFreeDir0).transpose() << "\n";
            std::cout << "dir from bound: " << dir_from_bound.transpose() << "\n";
            std::cout << "dir from free:  " << dir_from_free.transpose() << "\n---\n";
            ++n;
        }
    }

    std::cout << "failed / total = " << n << " / " << N << "\n";

}
