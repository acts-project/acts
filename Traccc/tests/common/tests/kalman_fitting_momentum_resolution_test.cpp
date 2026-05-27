/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <cmath>

// Local include(s).
#include "kalman_fitting_momentum_resolution_test.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// ROOT include(s).
#ifdef TRACCC_HAVE_ROOT
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#endif  // TRACCC_HAVE_ROOT

// System include(s).
#include <iostream>
#include <memory>
#include <stdexcept>

namespace traccc {

void KalmanFittingMomentumResolutionTests::consistency_tests(
    const edm::track_collection<default_algebra>::host::const_proxy_type& track,
    const edm::track_state_collection<default_algebra>::host&) const {

    // The nubmer of track states is supposed be equal to the number
    // of planes
    ASSERT_EQ(track.constituent_links().size(), std::get<11>(GetParam()));
}

void KalmanFittingMomentumResolutionTests::momentum_resolution_tests(
    [[maybe_unused]] std::string_view file_name) const {

#ifdef TRACCC_HAVE_ROOT

    // Open the file with the histograms.
    std::unique_ptr<TFile> ifile(TFile::Open(file_name.data(), "READ"));
    if ((!ifile) || ifile->IsZombie()) {
        throw std::runtime_error(std::string("Could not open file \"") +
                                 file_name.data() + "\"");
    }

    // Access the histogram.
    TH1* residual_qopT_hist = dynamic_cast<TH1*>(ifile->Get("res_qopT"));
    if (!residual_qopT_hist) {
        throw std::runtime_error(
            std::string("Could not access histogram residual_qopT in file \"") +
            file_name.data() + "\"");
    }

    // Function used for the fit.
    TF1 gaus{"gaus", "gaus", -5, 5};
    double fit_par[3];

    // Set the mean seed to 0
    gaus.SetParameters(1, 0.);
    gaus.SetParLimits(1, -1., 1.);

    // Set the standard deviation seed to 0.01
    gaus.SetParameters(2, 0.01f);
    gaus.SetParLimits(2, 0.f, 0.1f);

    auto res = residual_qopT_hist->Fit("gaus", "Q0S");

    gaus.GetParameters(&fit_par[0]);

    // Mean check
    EXPECT_NEAR(fit_par[1], 0, 0.05) << " Residual qopT mean value error";

    // Expected momentum resolution
    const std::size_t N = std::get<11>(GetParam());
    const auto smearing = std::get<15>(GetParam());
    const scalar epsilon = smearing[0u];
    const scalar Bz = std::get<13>(GetParam())[2u];
    const scalar spacing = std::get<12>(GetParam());
    const scalar L = static_cast<scalar>(N - 1u) * spacing;

    // Curvature (1/helix_radius) resolution from detector noise Eq. (35.61) of
    // PDG 2024
    const scalar dkres =
        epsilon / (L * L) * math::sqrt(720.f / static_cast<scalar>(N + 4));

    // Curvature (1/helix_radius) resolution from multiple scattering Eq.
    // (35.63) of PDG 2024
    // @NOTE: The calculation of the multiple scattering term is in work in
    // progress. Currently we are validating the momentum resolution without any
    // material on the modules, therefore, the multiple scattering contribution
    // to the momentum resolution is set to zero.
    const scalar dkms = 0.f;

    // Eq. (35.60) of PDG 2024
    const scalar dk = math::sqrt(dkres * dkres + dkms * dkms);

    // σ(q/pT) = σ(k/B) = σ(k)/B
    const scalar expected_sigma_qopT = dk / Bz;

    // Check if the standard deviation of the qopT residuals is within the
    // theoretically expected range.
    EXPECT_NEAR(fit_par[2], expected_sigma_qopT, expected_sigma_qopT * 0.025f)
        << " Residual qopT sigma value error";

#else
    std::cout << "Momentum resolution tests not performed without ROOT"
              << std::endl;
#endif  // TRACCC_HAVE_ROOT

    return;
}

void KalmanFittingMomentumResolutionTests::SetUp() {

    vecmem::host_memory_resource host_mr;

    const scalar offset = std::get<10>(GetParam());
    const unsigned int n_planes = std::get<11>(GetParam());
    const scalar spacing = std::get<12>(GetParam());

    std::vector<scalar> plane_positions;
    for (unsigned int i = 0; i < n_planes; i++) {
        plane_positions.push_back(offset * unit<scalar>::mm +
                                  static_cast<scalar>(i) * spacing *
                                      unit<scalar>::mm);
    }

    /// Plane alignment direction (aligned to x-axis)
    const vector3 bfield = std::get<13>(GetParam());

    const auto p = std::get<3>(GetParam());
    const detray::detail::helix<traccc::default_algebra> hlx(
        {0, 0, 0}, 0, {1, 0, 0}, -1.f / p, bfield);

    constexpr scalar rect_half_length = 500.f;
    constexpr detray::mask<detray::rectangle2D, traccc::default_algebra>
        rectangle{0u, rect_half_length, rect_half_length};

    detray::tel_det_config<default_algebra, detray::rectangle2D,
                           detray::detail::helix>
        tel_cfg(rectangle, hlx);
    tel_cfg.positions(plane_positions);
    tel_cfg.module_material(std::get<14>(GetParam()));
    tel_cfg.mat_thickness(thickness);
    tel_cfg.envelope(rect_half_length);

    // Create telescope detector
    auto [det, name_map] = build_telescope_detector(host_mr, tel_cfg);

    // Write detector file
    auto writer_cfg = detray::io::detector_writer_config{}
                          .format(detray::io::format::json)
                          .replace_files(true)
                          .write_material(true)
                          .path(std::get<0>(GetParam()));
    detray::io::write_detector(det, name_map, writer_cfg);
}

}  // namespace traccc
