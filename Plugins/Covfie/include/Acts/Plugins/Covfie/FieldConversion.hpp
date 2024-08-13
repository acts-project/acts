// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Covfie include(s)
#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/field_view.hpp>
#include <covfie/core/parameter_pack.hpp>

// Acts include(s)
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

namespace Acts::CovfiePlugin {

template <typename scalar_t, template <typename> class storage_t>
struct Converter{
    public:

    using Scalar = scalar_t;

    using Vector = typename covfie::vector::vector_d<Scalar, 3UL>;

    using Storage = storage_t<Vector>;

    using BuilderBackend = 
            covfie::backend::strided<covfie::vector::size3, Storage>;

    using InterpolatedField = covfie::field<covfie::backend::clamp<
        covfie::backend::affine<covfie::backend::linear<BuilderBackend>>>>;

    using ConstantField = covfie::field<
        covfie::backend::constant<Vector, Vector>>;

    /// @brief Creates a covfie field from an interpolated magnetic field.
    /// @param magneticField The acts interpolated magnetic field.
    /// @return An affine linear strided covfie field.
    InterpolatedField covfieField(
        const Acts::InterpolatedMagneticField& magneticField) const {
    Acts::MagneticFieldContext ctx;
    auto cache = magneticField.makeCache(ctx);
    return covfieFieldLinear(magneticField, cache, magneticField.getNBins(),
                            magneticField.getMin(), magneticField.getMax());
    }

    /// @brief Creates a covfie field from a constant B field.
    /// @param magneticField The acts constant magnetic field.
    /// @return A constant covfie field.
    ConstantField covfieField(const Acts::ConstantBField& magneticField) const {
        auto B = magneticField.getField();
        ConstantField field(
            covfie::make_parameter_pack(typename ConstantField::backend_t::configuration_t{
                static_cast<float>(B[0]), static_cast<float>(B[1]),
                static_cast<float>(B[2])}));
        return field;
    }

    /// @brief Creates a covfie field from a magnetic field provider by sampling it.
    /// The field must be defined within min (inclusive) and max (inclusive).
    /// @param magneticField The acts magnetic field provider.
    /// @param cache The acts cache.
    /// @param nBins 3D array of containing the number of bins for each axis.
    /// @param min (min_x, min_y, min_z)
    /// @param max (max_x, max_y, max_z)
    /// @return An affine linear strided covfie field.
    InterpolatedField covfieField(const Acts::MagneticFieldProvider& magneticField,
                                Acts::MagneticFieldProvider::Cache& cache,
                                const std::vector<std::size_t>& nBins,
                                const std::vector<Acts::ActsScalar>& min,
                                const std::vector<Acts::ActsScalar>& max) const {
    return covfieFieldLinear(magneticField, cache, nBins, min, max);
    }
    
    
    private:
    /// @brief Get the value of the interpolated field at a specific position in
    /// min (inclusive) max (inclusive).
    /// @note The unchecked lookup allows the access to the field value at the max edges
    /// since the domain is typically min (inclusive) max (exclusive).
    /// @param magneticField the interpolated magnetic field.
    /// @param cache the magnetic field cache.
    /// @param position the position of the field to look up.
    /// @return the field value at the given position.
    Acts::Vector3 getFieldEdgeInclusive(
        const Acts::InterpolatedMagneticField& magneticField,
        [[maybe_unused]] Acts::InterpolatedMagneticField::Cache& cache,
        const Acts::Vector3& position) const {
    // Check if position in within [min; max].
    bool inBounds = position[0] <= magneticField.getMax()[0] &&
                    position[1] <= magneticField.getMax()[1] &&
                    position[2] <= magneticField.getMax()[2] &&
                    position[0] >= magneticField.getMin()[0] &&
                    position[1] >= magneticField.getMin()[1] &&
                    position[2] >= magneticField.getMin()[2];
    if (!inBounds) {
        throw std::runtime_error{
            "Field lookup failure (position not in [min, max])"};
    }
    return magneticField.getFieldUnchecked(position);
    }

    /// @brief Get the value of the field at a specific position of a general magnetic field.
    /// @param magneticField the interpolated magnetic field.
    /// @param cache the magnetic field cache.
    /// @param position the position of the field to look up.
    /// @return the field value at the given position.
    Acts::Vector3 getFieldEdgeInclusive(
        const Acts::MagneticFieldProvider& magneticField,
        Acts::MagneticFieldProvider::Cache& cache, const Acts::Vector3& position) const {
    auto lookupResult = magneticField.getField(position, cache);
    if (!lookupResult.ok()) {
        throw std::runtime_error{"Field lookup failure"};
    }
    return *lookupResult;
    }

    /// @brief Creates a strided covfie field that stores the values of the magnetic field in the volume given by min and max using a fixed sample spacing (determined by nBins).
    /// @param magneticField The acts magnetic field.
    /// @param cache The acts cache.
    /// @param nBins 3D array of containing the number of bins for each axis.
    /// @param min (min_x, min_y, min_z)
    /// @param max (max_x, max_y, max_z)
    /// @return A strided covfie field.
    template <typename magnetic_field_t>
    auto newBuilder(const magnetic_field_t& magneticField,
                    typename magnetic_field_t::Cache& cache,
                    const std::vector<std::size_t>& nBins,
                    const std::vector<Acts::ActsScalar>& min,
                    const std::vector<Acts::ActsScalar>& max) const {
    using Field = covfie::field<BuilderBackend>;

    Field field(covfie::make_parameter_pack(
        typename Field::backend_t::configuration_t{nBins[0], nBins[1], nBins[2]}));

    typename Field::view_t view(field);

    std::array<double, 3> sampleSpacing = {(max[0] - min[0]) / (nBins[0] - 1),
                                            (max[1] - min[1]) / (nBins[1] - 1),
                                            (max[2] - min[2]) / (nBins[2] - 1)};

    for (std::size_t x = 0; x < nBins[0]; x++) {
        for (std::size_t y = 0; y < nBins[1]; y++) {
        for (std::size_t z = 0; z < nBins[2]; z++) {
            auto position = Acts::Vector3{x * sampleSpacing[0] + min[0],
                                        y * sampleSpacing[1] + min[1],
                                        z * sampleSpacing[2] + min[2]};

            auto result = getFieldEdgeInclusive(magneticField, cache, position);

            typename Field::view_t::output_t& p = view.at(x, y, z);
            p[0] = static_cast<float>(result[0]);
            p[1] = static_cast<float>(result[1]);
            p[2] = static_cast<float>(result[2]);
        }
        }
    }
    return field;
    }

    /// @brief Generate the affine covfie configuration (scaling and rotation) given the size of the field (min and max)
    /// @param nBins 3D array of containing the number of bins for each axis.
    /// @param min (min_x, min_y, min_z)
    /// @param max (max_x, max_y, max_z)
    /// @return The affine field configuration.
    template <typename backend_t>
    typename backend_t::configuration_t affineConfiguration(
        const std::vector<std::size_t>& nBins,
        const std::vector<Acts::ActsScalar>& min,
        const std::vector<Acts::ActsScalar>& max) const {
    auto scaling = covfie::algebra::affine<3>::scaling(
        static_cast<float>((nBins[0] - 1) / (max[0] - min[0])),
        static_cast<float>((nBins[1] - 1) / (max[1] - min[1])),
        static_cast<float>((nBins[2] - 1) / (max[2] - min[2])));

    auto translation = covfie::algebra::affine<3>::translation(
        static_cast<float>(-min[0]), static_cast<float>(-min[1]),
        static_cast<float>(-min[2]));

    return {scaling * translation};
    }

    /// @brief Uses std::nextafter to generates a clamp backend
    /// configuration where arguments min and max hold floating point values.
    /// @param min (min_x, min_y, min_z)
    /// @param max (max_x, max_y, max_z)
    /// @return The clamp field configuration.
    template <typename backend_t>
    typename backend_t::configuration_t clampConfigurationFloat(
        const std::vector<Acts::ActsScalar>& min,
        const std::vector<Acts::ActsScalar>& max) const {
    return {{std::nextafter(static_cast<float>(min[0]),
                            std::numeric_limits<float>::infinity()),
            std::nextafter(static_cast<float>(min[1]),
                            std::numeric_limits<float>::infinity()),
            std::nextafter(static_cast<float>(min[2]),
                            std::numeric_limits<float>::infinity())},
            {std::nextafter(static_cast<float>(max[0]),
                            -std::numeric_limits<float>::infinity()),
            std::nextafter(static_cast<float>(max[1]),
                            -std::numeric_limits<float>::infinity()),
            std::nextafter(static_cast<float>(max[2]),
                            -std::numeric_limits<float>::infinity())}};
    }

    /// @brief Creates a covfie field from a generic magnetic field.
    /// @param magneticField The generic magnetic field.
    /// @param cache The cache.
    /// @param nBins 3D array of containing the number of bins for each axis.
    /// @param min (min_x, min_y, min_z)
    /// @param max (max_x, max_y, max_z)
    /// @return A clamp affine linear strided covfie field.
    template <typename magnetic_field_t>
    InterpolatedField covfieFieldLinear(const magnetic_field_t& magneticField,
                                        typename magnetic_field_t::Cache& cache,
                                        const std::vector<std::size_t>& nBins,
                                        const std::vector<Acts::ActsScalar>& min,
                                        const std::vector<Acts::ActsScalar>& max) const {
    auto builder = newBuilder(magneticField, cache, nBins, min, max);
    InterpolatedField field(covfie::make_parameter_pack(
        clampConfigurationFloat<typename InterpolatedField::backend_t>(min, max),
        affineConfiguration<typename InterpolatedField::backend_t::backend_t>(nBins, min,
                                                                    max),
        typename InterpolatedField::backend_t::backend_t::backend_t::configuration_t{},
        builder.backend()));

    return field;
    }
};

}  // namespace Acts::CovfiePlugin
