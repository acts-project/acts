/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <optional>
#include <string>
#include <traccc/definitions/primitives.hpp>
#include <traccc/edm/particle.hpp>
#include <vector>

namespace traccc {
/**
 * @brief Abstract track matching descriptor.
 *
 * This abstract class represents a methodology for matching seeds to
 * particles.
 */
struct track_matcher {
    /**
     * @brief Return a human-readable name for this matching method.
     */
    virtual std::string get_name() const = 0;

    /**
     * @brief Determine the most likely particle, if any, for a given seed.
     *
     * @param v A vector of sets of particle identifiers, such that the outer
     * vector represents the list of spacepoints in a seed, and the inner sets
     * represent the particles that match that spacepoint (usually only one).
     *
     * @return The matched particle identifier, or nothing if no particle
     * matches.
     */
    virtual std::optional<uint64_t> operator()(
        const std::vector<std::vector<uint64_t>>& v) const = 0;
    /**
     * @brief Default destructor
     */
    virtual ~track_matcher() = default;
};

/**
 * @brief Matcher based on minimum majority size.
 *
 * This matcher matches a seed to the track with the greatest number of
 * matching spacepoints, with a minimum percentage-defined size of the
 * majority.
 */
struct stepped_percentage : track_matcher {
    /**
     * @brief Construct a new particle matcher.
     *
     * @param r The minimum ratio in [0, 1] of spacepoints that must belong to
     * the same particle.
     */
    stepped_percentage(scalar r);

    virtual std::string get_name() const override final;

    virtual std::optional<uint64_t> operator()(
        const std::vector<std::vector<uint64_t>>&) const override final;

    private:
    scalar m_min_ratio;
};

/**
 * @brief Matcher requiring exact particle matching.
 *
 * This matcher matches a seed to a track iff all spacepoints in the seed map
 * onto the same particle.
 */
struct exact : track_matcher {
    exact();

    virtual std::string get_name() const override final;

    virtual std::optional<uint64_t> operator()(
        const std::vector<std::vector<uint64_t>>&) const override final;
};
}  // namespace traccc
