#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/DiskSurface.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/EventData/detail/TransformationBoundToFree.hpp"

#include "Acts/Utilities/detail/gaussian_mixture_helpers.hpp"

int main() {
    std::srand((unsigned int) std::time(0));

    //auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(Acts::Vector3{0,0,0}, Acts::Vector3{1,0,0});

    auto disk_bounds = Acts::DiscBounds(

    for(int i=0; i<100'000; ++i) {

        Acts::FreeVector pars1 = Acts::FreeVector::Zero();
        Acts::FreeVector pars2 = Acts::FreeVector::Zero();

        pars1.segment<3>(Acts::eFreeDir0) = Acts::Vector3{1,0,0}; //Acts::Vector3::Random().normalized();
        pars2.segment<3>(Acts::eFreeDir0) = Acts::Vector3{0,1,0}; //Acts::Vector3::Random().normalized();

        const auto bound1 = *Acts::detail::transformFreeToBoundParameters(pars1, *surface, Acts::GeometryContext{});
        const auto bound2 = *Acts::detail::transformFreeToBoundParameters(pars2, *surface, Acts::GeometryContext{});

        std::array<std::tuple<double, Acts::BoundVector, std::optional<Acts::BoundSymMatrix>>, 2> a;

        a[0] = {0.5, bound1, Acts::BoundSymMatrix::Identity().eval() };
        a[1] = {0.5, bound2, Acts::BoundSymMatrix::Identity().eval() };

        const auto [res_bound, cov] = Acts::detail::combineBoundGaussianMixture(a.begin(), a.end());

        const auto res_free = Acts::detail::transformBoundToFreeParameters(*surface, Acts::GeometryContext{}, res_bound);

        const auto dir_before = (0.5 * (pars1.segment<3>(Acts::eFreeDir0) + pars2.segment<3>(Acts::eFreeDir0))).normalized();
        const auto dir_after = res_free.segment<3>(Acts::eFreeDir0).normalized();

        const double res = (dir_before - dir_after).array().abs().sum();

        if( res > 1.e-4 ) {
            std::cout << "res: " << res << "\n";
            std::cout << "first before:  " << pars1.segment<3>(Acts::eFreeDir0).transpose() << "\n";
            std::cout << "second before: " << pars2.segment<3>(Acts::eFreeDir0).transpose() << "\n";
            std::cout << "combinde before: " << dir_before.transpose() << "\n";
            std::cout << "combinde after: "  << dir_after.transpose() << "\n";
        }
    }


}
