// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "traccc/edm/measurement.hpp"
#include "traccc/edm/spacepoint.hpp"
#include "traccc/edm/seed.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <type_traits>

namespace ActsExamples::Traccc::Common::Util {

template <typename S>
struct Hasher{
    std::size_t operator()(const S& s) const noexcept {
        return hash(s);
    }
    private:
    std::hash<S> hash;
};

template <>
struct Hasher<traccc::measurement>{
    std::size_t operator()(const traccc::measurement& s) const noexcept {
        return s.measurement_id;
    }
};

template <>
struct Hasher<traccc::spacepoint>{
    std::size_t operator()(const traccc::spacepoint& s) const noexcept {
        return Hasher<traccc::measurement>{}(s.meas);
    }
};

template <>
struct Hasher<traccc::seed>{
    std::size_t operator()(const traccc::seed& s) const noexcept {
        return s.spB_link ^ s.spM_link ^ s.spT_link;
    }
};

template <typename T>
struct Equals{
    std::size_t operator()(const T& lhs, const T& rhs) const {
        return eq(lhs, rhs);
    }
    private:
    std::equal_to<T> eq;
};

template <>
struct Equals<traccc::seed>{
    std::size_t operator()(const traccc::seed& lhs, const traccc::seed& rhs) const {
        return lhs.spB_link == rhs.spB_link &&
        lhs.spM_link == rhs.spM_link &&
        lhs.spT_link == rhs.spT_link &&
        lhs.weight == rhs.weight &&
        lhs.z_vertex == rhs.z_vertex;
    }
};

}