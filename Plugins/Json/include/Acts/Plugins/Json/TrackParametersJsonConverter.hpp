#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/SurfaceJsonConverter.hpp"

#include <nlohmann/json.hpp>

template <class parameters_t>
concept TrackParameters = 
    Acts::FreeTrackParametersConcept<parameters_t> || 
    Acts::BoundTrackParametersConcept<parameters_t>;

template <class parameters_t>
concept IsGenericBound = 
    std::same_as<
        parameters_t, 
        Acts::GenericBoundTrackParameters<
            typename parameters_t::ParticleHypothesis>>;

namespace nlohmann {
    
    template <TrackParameters parameters_t>
    struct adl_serializer<parameters_t> {
        using CovarianceMatrix = typename parameters_t::CovarianceMatrix;
        static void to_json(nlohmann::json& j, const parameters_t& t) {
            j["direction"] = t.direction();
            j["qOverP"] = t.qOverP();
            j["particleHypothesis"] = t.particleHypothesis().absolutePdg();
            
            j["covariance"];
            if (t.covariance().has_value()) {
                auto cov = t.covariance().value();
                constexpr unsigned int size = cov.rows();
                std::array<Acts::ActsScalar, size*size> covData;
                for (size_t n = 0; n < size; ++n) {
                    for (size_t m = 0; m < size; ++m) {
                        covData[n * size + m] = cov(n, m);
                    }
                }
                j["covariance"] = covData;
            }
            if constexpr (IsGenericBound<parameters_t>) {
                Acts::GeometryContext gctx;
                j["position"] = t.fourPosition(gctx);

                j["referenceSurface"] = Acts::SurfaceJsonConverter::toJson(gctx, t.referenceSurface());
            }
            else {
                j["position"] = t.fourPosition();
            }
        }
    
        static parameters_t from_json(const nlohmann::json& j) {
            std::array<Acts::ActsScalar, 4> posData = j.at("position");
            Acts::Vector4 position(posData[0], posData[1], posData[2], posData[3]);

            std::array<Acts::ActsScalar, 3> dirData = j.at("direction");
            Acts::Vector3 direction(dirData[0], dirData[1], dirData[2]);

            Acts::ActsScalar qOverP = j.at("qOverP");
            Acts::PdgParticle absPdg = j.at("particleHypothesis");
            
            std::optional<CovarianceMatrix> cov;
            if (j.at("covariance").is_null()) {
                cov = std::nullopt;
            }
            else {
                CovarianceMatrix mat;
                constexpr unsigned int size = mat.rows();
                std::array<Acts::ActsScalar, size*size> covData = j.at("covariance");
                for (size_t n = 0; n < size; ++n) {
                    for (size_t m = 0; m < size; ++m) {
                        mat(n, m) = covData[n * size + m];
                    }
                }
                cov.emplace(std::move(mat));
            }

            typename parameters_t::ParticleHypothesis particle(absPdg);

            if constexpr (IsGenericBound<parameters_t>) {
                Acts::GeometryContext gctx;
                auto referenceSurface = Acts::SurfaceJsonConverter::fromJson(j.at("referenceSurface"));

                auto res = parameters_t::create(
                    referenceSurface, gctx, position, direction, qOverP, cov, particle);

                if (!res.ok()) {
                    throw std::invalid_argument("Invalid bound track parameters");
                }
                return res.value();
            }
            else {
                return parameters_t(position, direction, qOverP, cov, particle);
            }
        }
    };

    template <TrackParameters parameters_t>
    struct adl_serializer<std::shared_ptr<parameters_t>> {
        using CovarianceMatrix = typename parameters_t::CovarianceMatrix;
        static void to_json(nlohmann::json& j, const parameters_t& t) {
            j = *t;
        }
    
        static std::shared_ptr<parameters_t> from_json(const nlohmann::json& j) {
            return std::make_shared<parameters_t>(j.get<parameters_t>());
        }
    };

}   // namespace nlohmann

namespace Acts {

    template <TrackParameters parameters_t>
    void to_json(nlohmann::json& j, const std::unique_ptr<parameters_t>& t) {
        if (t == nullptr) {
            return;
        }
        j = *t;
    }

    template <TrackParameters parameters_t>
    void to_json(nlohmann::json& j, const std::unique_ptr<const parameters_t>& t) {
        if (t == nullptr) {
            return;
        }
        j = *t;
    }

    template <TrackParameters parameters_t>
    void to_json(nlohmann::json& j, const std::shared_ptr<parameters_t>& t) {
        if (t == nullptr) {
            return;
        }
        j = *t;
    }

    template <TrackParameters parameters_t>
    void to_json(nlohmann::json& j, const std::shared_ptr<const parameters_t>& t) {
        if (t == nullptr) {
            return;
        }
        j = *t;
    }

}  // namespace Acts
