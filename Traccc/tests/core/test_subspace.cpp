/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "traccc/utils/subspace.hpp"

using namespace traccc;

TEST(subspace, init) {
    const std::array<detray::dindex_type<default_algebra>, 2> indices{1, 2};
    subspace<default_algebra, 6, 2> s{indices};

    ASSERT_EQ(s.get_index(0), indices[0]);
    ASSERT_EQ(s.get_index(1), indices[1]);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    const auto H = s.projector<2>();

    ASSERT_EQ(getter::element(H, 0, 0), 0);
    ASSERT_EQ(getter::element(H, 0, 1), 1);
    ASSERT_EQ(getter::element(H, 0, 2), 0);
    ASSERT_EQ(getter::element(H, 0, 3), 0);
    ASSERT_EQ(getter::element(H, 0, 4), 0);
    ASSERT_EQ(getter::element(H, 0, 5), 0);
    ASSERT_EQ(getter::element(H, 1, 0), 0);
    ASSERT_EQ(getter::element(H, 1, 1), 0);
    ASSERT_EQ(getter::element(H, 1, 2), 1);
    ASSERT_EQ(getter::element(H, 1, 3), 0);
    ASSERT_EQ(getter::element(H, 1, 4), 0);
    ASSERT_EQ(getter::element(H, 1, 5), 0);
}

TEST(subspace, set_index) {
    const std::array<detray::dindex_type<default_algebra>, 2> indices{1, 2};
    subspace<default_algebra, 6, 2> s{indices};

    s.set_index(0, 5);

    ASSERT_EQ(s.get_index(0), 5);
    ASSERT_EQ(s.get_index(1), indices[1]);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    s.set_index(1, 4);

    ASSERT_EQ(s.get_index(0), 5);
    ASSERT_EQ(s.get_index(1), 4);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    s.set_index(0, 3);

    ASSERT_EQ(s.get_index(0), 3);
    ASSERT_EQ(s.get_index(1), 4);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    s.set_index(1, 0);

    ASSERT_EQ(s.get_index(0), 3);
    ASSERT_EQ(s.get_index(1), 0);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    const auto H = s.projector<2>();

    ASSERT_EQ(getter::element(H, 0, 0), 0);
    ASSERT_EQ(getter::element(H, 0, 1), 0);
    ASSERT_EQ(getter::element(H, 0, 2), 0);
    ASSERT_EQ(getter::element(H, 0, 3), 1);
    ASSERT_EQ(getter::element(H, 0, 4), 0);
    ASSERT_EQ(getter::element(H, 0, 5), 0);
    ASSERT_EQ(getter::element(H, 1, 0), 1);
    ASSERT_EQ(getter::element(H, 1, 1), 0);
    ASSERT_EQ(getter::element(H, 1, 2), 0);
    ASSERT_EQ(getter::element(H, 1, 3), 0);
    ASSERT_EQ(getter::element(H, 1, 4), 0);
    ASSERT_EQ(getter::element(H, 1, 5), 0);
}

TEST(subspace, set_sign) {
    const std::array<detray::dindex_type<default_algebra>, 2> indices{1, 2};
    subspace<default_algebra, 6, 2> s{indices};

    s.set_sign(0, true);

    ASSERT_EQ(s.get_index(0), 1);
    ASSERT_EQ(s.get_index(1), 2);
    ASSERT_EQ(s.get_sign(0), true);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    s.set_sign(1, true);

    ASSERT_EQ(s.get_index(0), 1);
    ASSERT_EQ(s.get_index(1), 2);
    ASSERT_EQ(s.get_sign(0), true);
    ASSERT_EQ(s.get_sign(1), true);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    const auto H1 = s.projector<2>();

    ASSERT_EQ(getter::element(H1, 0, 0), 0);
    ASSERT_EQ(getter::element(H1, 0, 1), -1);
    ASSERT_EQ(getter::element(H1, 0, 2), 0);
    ASSERT_EQ(getter::element(H1, 0, 3), 0);
    ASSERT_EQ(getter::element(H1, 0, 4), 0);
    ASSERT_EQ(getter::element(H1, 0, 5), 0);
    ASSERT_EQ(getter::element(H1, 1, 0), 0);
    ASSERT_EQ(getter::element(H1, 1, 1), 0);
    ASSERT_EQ(getter::element(H1, 1, 2), -1);
    ASSERT_EQ(getter::element(H1, 1, 3), 0);
    ASSERT_EQ(getter::element(H1, 1, 4), 0);
    ASSERT_EQ(getter::element(H1, 1, 5), 0);

    s.set_sign(0, false);

    ASSERT_EQ(s.get_index(0), 1);
    ASSERT_EQ(s.get_index(1), 2);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), true);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    const auto H2 = s.projector<2>();

    ASSERT_EQ(getter::element(H2, 0, 0), 0);
    ASSERT_EQ(getter::element(H2, 0, 1), 1);
    ASSERT_EQ(getter::element(H2, 0, 2), 0);
    ASSERT_EQ(getter::element(H2, 0, 3), 0);
    ASSERT_EQ(getter::element(H2, 0, 4), 0);
    ASSERT_EQ(getter::element(H2, 0, 5), 0);
    ASSERT_EQ(getter::element(H2, 1, 0), 0);
    ASSERT_EQ(getter::element(H2, 1, 1), 0);
    ASSERT_EQ(getter::element(H2, 1, 2), -1);
    ASSERT_EQ(getter::element(H2, 1, 3), 0);
    ASSERT_EQ(getter::element(H2, 1, 4), 0);
    ASSERT_EQ(getter::element(H2, 1, 5), 0);
}

TEST(subspace, set_invalid) {
    const std::array<detray::dindex_type<default_algebra>, 2> indices{1, 2};
    subspace<default_algebra, 6, 2> s{indices};

    s.set_invalid(1);

    ASSERT_EQ(s.get_index(0), 1);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), false);

    const auto H1 = s.projector<2>();

    ASSERT_EQ(getter::element(H1, 0, 0), 0);
    ASSERT_EQ(getter::element(H1, 0, 1), 1);
    ASSERT_EQ(getter::element(H1, 0, 2), 0);
    ASSERT_EQ(getter::element(H1, 0, 3), 0);
    ASSERT_EQ(getter::element(H1, 0, 4), 0);
    ASSERT_EQ(getter::element(H1, 0, 5), 0);
    ASSERT_EQ(getter::element(H1, 1, 0), 0);
    ASSERT_EQ(getter::element(H1, 1, 1), 0);
    ASSERT_EQ(getter::element(H1, 1, 2), 0);
    ASSERT_EQ(getter::element(H1, 1, 3), 0);
    ASSERT_EQ(getter::element(H1, 1, 4), 0);
    ASSERT_EQ(getter::element(H1, 1, 5), 0);

    s.set_index(1, 5);

    ASSERT_EQ(s.get_index(0), 1);
    ASSERT_EQ(s.get_index(1), 5);
    ASSERT_EQ(s.get_sign(0), false);
    ASSERT_EQ(s.get_sign(1), false);
    ASSERT_EQ(s.get_valid(0), true);
    ASSERT_EQ(s.get_valid(1), true);

    const auto H2 = s.projector<2>();

    ASSERT_EQ(getter::element(H2, 0, 0), 0);
    ASSERT_EQ(getter::element(H2, 0, 1), 1);
    ASSERT_EQ(getter::element(H2, 0, 2), 0);
    ASSERT_EQ(getter::element(H2, 0, 3), 0);
    ASSERT_EQ(getter::element(H2, 0, 4), 0);
    ASSERT_EQ(getter::element(H2, 0, 5), 0);
    ASSERT_EQ(getter::element(H2, 1, 0), 0);
    ASSERT_EQ(getter::element(H2, 1, 1), 0);
    ASSERT_EQ(getter::element(H2, 1, 2), 0);
    ASSERT_EQ(getter::element(H2, 1, 3), 0);
    ASSERT_EQ(getter::element(H2, 1, 4), 0);
    ASSERT_EQ(getter::element(H2, 1, 5), 1);
}
