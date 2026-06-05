/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/options/track_seeding.hpp"

#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <format>

namespace traccc::opts {

/// Convenience namespace shorthand
namespace po = boost::program_options;

track_seeding::track_seeding() : interface("Track Seeding Options") {

    m_desc.add_options()(
        "seedfinder-z-range",
        po::value(&m_z_range)->value_name("MIN:MAX")->default_value(m_z_range),
        "Spacepoint Z range [mm]");
    m_desc.add_options()(
        "seedfinder-r-range",
        po::value(&m_r_range)->value_name("MIN:MAX")->default_value(m_r_range),
        "Spacepoint R range [mm]");
    m_desc.add_options()("seedfinder-vertex-range",
                         po::value(&m_vertex_range)
                             ->value_name("MIN:MAX")
                             ->default_value(m_vertex_range),
                         "Vertex Z range [mm]");

    m_desc.add_options()(
        "seedfinder-minPt",
        po::value(&m_seedfinder.minPt)
            ->default_value(m_seedfinder.minPt / unit<float>::GeV),
        "Minimum track momentum [GeV]");

    m_desc.add_options()("seedfinder-cotThetaMax",
                         po::value(&m_seedfinder.cotThetaMax)
                             ->default_value(m_seedfinder.cotThetaMax),
                         "Maximum cotangent of theta angle [unitless]");
    m_desc.add_options()(
        "seedfinder-deltaR-range",
        po::value(&m_delta_r_range)->default_value(m_delta_r_range),
        "Radial distance between measurements [mm]");

    m_desc.add_options()(
        "seedfinder-impactMax",
        po::value(&m_seedfinder.impactMax)
            ->default_value(m_seedfinder.impactMax / unit<float>::mm),
        "Maximum impact parameter [mm]");

    m_desc.add_options()(
        "seedfinder-deltaZMax",
        po::value(&m_seedfinder.deltaZMax)
            ->default_value(m_seedfinder.deltaZMax / unit<float>::mm),
        "Maximum Z distance between measurements [mm]");

    m_desc.add_options()("seedfinder-sigmaScattering",
                         po::value(&m_seedfinder.sigmaScattering)
                             ->default_value(m_seedfinder.sigmaScattering),
                         "Scattering angle sigma [unitless]");
    m_desc.add_options()(
        "seedfinder-maxPtScattering",
        po::value(&m_seedfinder.maxPtScattering)
            ->default_value(m_seedfinder.maxPtScattering / unit<float>::GeV),
        "Upper pT limit for scattering [GeV]");

    m_desc.add_options()("seedfinder-maxSeedsPerSpM",
                         po::value(&m_seedfinder.maxSeedsPerSpM)
                             ->default_value(m_seedfinder.maxSeedsPerSpM),
                         "Maximum number of seeds per middle space point");
    m_desc.add_options()(
        "seedfinder-bFieldInZ",
        po::value(&m_seedfinder.bFieldInZ)
            ->default_value(m_seedfinder.bFieldInZ / unit<float>::T),
        "B-field in Z direction [T]");
    m_desc.add_options()(
        "seedfinder-deltaInvHelixDiameter",
        po::value(&m_seedfilter.deltaInvHelixDiameter)
            ->default_value(m_seedfilter.deltaInvHelixDiameter *
                            unit<float>::mm),
        "Inverted radius delta for compatible seeds [mm^-1]");
    m_desc.add_options()("seedfinder-impactWeightFactor",
                         po::value(&m_seedfilter.impactWeightFactor)
                             ->default_value(m_seedfilter.impactWeightFactor),
                         "Weight factor for impact parameter [unitless]");
    m_desc.add_options()("seedfinder-compatSeedWeight",
                         po::value(&m_seedfilter.compatSeedWeight)
                             ->default_value(m_seedfilter.compatSeedWeight),
                         "Weight per compatible seed [unitless]");
    m_desc.add_options()("seedfinder-compatSeedLimit",
                         po::value(&m_seedfilter.compatSeedLimit)
                             ->default_value(m_seedfilter.compatSeedLimit),
                         "Maximum weighted compatible seeds [cardinal]");
}

track_seeding::operator seedfinder_config() const {

    return m_seedfinder;
}

track_seeding::operator seedfilter_config() const {

    return m_seedfilter;
}

track_seeding::operator spacepoint_grid_config() const {

    return {m_seedfinder};
}

track_seeding::operator vector3() const {

    return {0.f, 0.f, m_seedfinder.bFieldInZ};
}

void track_seeding::read(const po::variables_map&) {

    m_seedfinder.zMin = m_z_range[0] * unit<float>::mm;
    m_seedfinder.zMax = m_z_range[1] * unit<float>::mm;
    m_seedfinder.rMin = m_r_range[0] * unit<float>::mm;
    m_seedfinder.rMax = m_r_range[1] * unit<float>::mm;
    m_seedfinder.collisionRegionMin = m_vertex_range[0] * unit<float>::mm;
    m_seedfinder.collisionRegionMax = m_vertex_range[1] * unit<float>::mm;
    m_seedfinder.deltaRMin = m_delta_r_range[0] * unit<float>::mm;
    m_seedfinder.deltaRMax = m_delta_r_range[1] * unit<float>::mm;

    m_seedfinder.minPt *= unit<float>::GeV;
    m_seedfinder.impactMax *= unit<float>::mm;
    m_seedfinder.maxPtScattering *= unit<float>::GeV;
    m_seedfinder.bFieldInZ *= unit<float>::T;

    m_seedfilter.deltaInvHelixDiameter /= unit<float>::mm;

    m_seedfinder.setup();
}

std::unique_ptr<configuration_printable> track_seeding::as_printable() const {

    auto cat = std::make_unique<configuration_category>(m_description);

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Spacepoint Z range",
        std::format("[{:.1f} - {:.1f}] mm", m_seedfinder.zMin / unit<float>::mm,
                    m_seedfinder.zMax / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Spacepoint R range",
        std::format("[{:.2f} - {:.2f}] mm", m_seedfinder.rMin / unit<float>::mm,
                    m_seedfinder.rMax / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Vertex Z range",
        std::format("[{:.2f} - {:.2f}] mm",
                    m_seedfinder.collisionRegionMin / unit<float>::mm,
                    m_seedfinder.collisionRegionMax / unit<float>::mm)));

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Minimum track momentum",
        std::format("{:.2f} GeV", m_seedfinder.minPt / unit<float>::GeV)));

    float theta = std::atan(1.f / m_seedfinder.cotThetaMax);
    float eta = -std::log(std::tan(theta / 2.f));

    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum cotangent of theta angle",
        std::format("{:.4f} (eta = {:.2f})", m_seedfinder.cotThetaMax, eta)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Radial distance between measurements",
        std::format("[{:.2f} - {:.2f}] mm",
                    m_seedfinder.deltaRMin / unit<float>::mm,
                    m_seedfinder.deltaRMax / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum impact parameter",
        std::format("{:.2f} mm", m_seedfinder.impactMax / unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Scattering angle sigma",
        std::format("{:.2f}", m_seedfinder.sigmaScattering)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Upper pT limit for scattering",
        std::format("{:.2f} GeV",
                    m_seedfinder.maxPtScattering / unit<float>::GeV)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum seeds per middle space point",
        std::to_string(m_seedfinder.maxSeedsPerSpM)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "B-field in Z direction",
        std::format("{:.2f} T", m_seedfinder.bFieldInZ / unit<float>::T)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Inverted radius delta for compatible seeds",
        std::format("{:.5f} mm^-1",
                    m_seedfilter.deltaInvHelixDiameter * unit<float>::mm)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Weight factor for impact parameter",
        std::format("{:.2f}", m_seedfilter.impactWeightFactor)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Weight per compatible seed",
        std::format("{:.2f}", m_seedfilter.compatSeedWeight)));
    cat->add_child(std::make_unique<configuration_kv_pair>(
        "Maximum weighted compatible seed",
        std::format("{:d}", m_seedfilter.compatSeedLimit)));

    return cat;
}
}  // namespace traccc::opts
