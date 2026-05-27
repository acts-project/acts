/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <string>
#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/particle.hpp>

namespace traccc {
/**
 * @brief Abstract track filtering descriptor.
 *
 * This abstract class represents a methodology for selecting particles that
 * are relevant to the current track reconstruction. Allows irrelevant tracks
 * to be omitted from efficiency calculations.
 */
struct track_filter {
    /**
     * @brief Return a human-readable name for this filtering method.
     */
    virtual std::string get_name() const = 0;

    /**
     * @brief Determine whether a given particle is relevant to tracking.
     *
     * @param p The particle in question.
     *
     * @return true if the particle passes the cut.
     * @return false if the particle does not pass the cut.
     */
    virtual bool operator()(const particle &p) const = 0;

    /**
     * @brief Default destructor
     */
    virtual ~track_filter() = default;
};

/**
 * @brief Simple η and pT cut for charged particles.
 *
 * Describes a cut requiring the absolute value of η to be smaller than a given
 * value, and for the pT value to be greater than a given threshold. Also
 * requires the particle to be charged.
 */
struct simple_charged_eta_pt_cut : track_filter {
    /**
     * @brief Construct a new cut object.
     *
     * @param eta The maximum value for |η|
     * @param pT The minimum value for pT
     */
    simple_charged_eta_pt_cut(scalar eta, scalar pT);

    virtual std::string get_name() const override final;

    virtual bool operator()(const particle &) const override final;

    private:
    scalar m_eta, m_pT;
};
}  // namespace traccc
