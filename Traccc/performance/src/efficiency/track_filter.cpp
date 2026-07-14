/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <string>
#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/particle.hpp>
#include <traccc/efficiency/track_filter.hpp>

namespace traccc {
simple_charged_eta_pt_cut::simple_charged_eta_pt_cut(scalar eta, scalar pt)
    : m_eta(eta), m_pT(pt) {}

std::string simple_charged_eta_pt_cut::get_name() const {
    char buffer[512];
    snprintf(buffer, 512, "Charged with |η| ≤ %.2f and pT ≥ %.3f GeV", m_eta,
             m_pT);
    return std::string(buffer);
}

bool simple_charged_eta_pt_cut::operator()(const particle& p) const {
    const scalar eta = vector::eta(p.momentum);
    const scalar pT = vector::perp(p.momentum);

    return p.charge != 0 && std::abs(eta) <= m_eta && pT >= m_pT;
}
}  // namespace traccc
