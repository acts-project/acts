// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"

#include <set>

namespace Acts{
    /** @brief Factory helper class to construct volume or surface bounds where the constructed bounds
     *         are cached inside the factory and if the same bound parameters are requested at a later
     *         stage the factory automatically returns the cached bounds. This provides a simple sharing
     *         mechanism of the same bounds across multiple surfaces / volumes
       */
    
    
    namespace detail {       
        /** @brief Concept to define the minimal requirements on the bounds to apply the deduplication mechanism
         *         using the factory. */
        template <typename BoundsType_t> 
            concept ComparableBoundConcept = requires(const BoundsType_t& bounds) {
            /** @brief Getter function to distinguish the various bound types (e.g box vs. cylinder) */
            {bounds.type()};
            /** @brief Getter function returning all defining parameters of the bounds as an std::vector */
            {bounds.values()} ->std::same_as<std::vector<double>>;
        };
    }
    
    /** @brief  */
    template<detail::ComparableBoundConcept BoundsType_t>
        class BoundFactory {
            public:
                /** @brief Empty default constructor */
                BoundFactory() = default;
                /** @brief Delete the copy constructor */
                BoundFactory(const BoundFactory& other) = delete;
                /** @brief Delete copy assignment */
                BoundFactory& operator=(const BoundFactory& other) = delete;
                /** @brief Pass externally constructed bounds to the factory and run the deduplication
                 *         mechanism on them 
                 *  @param bounds: Pointer to the bounds to deduplicated */
                template <typename BoundsImpl_t>
                    std::shared_ptr<BoundsImpl_t> insert(const std::shared_ptr<BoundsImpl_t>& bounds)
                        requires(std::is_base_of_v<BoundsType_t, BoundsImpl_t>) {
                        assert(bounds);
                        return std::dynamic_pointer_cast<BoundsImpl_t>(*m_boundSet.insert(bounds).first);

                    } 
                template<typename BoundsImpl_t,
                         class... argList> 
                    std::shared_ptr<const BoundsImpl_t> makeBounds(argList... args)
                        requires(std::is_base_of_v<BoundsType_t, BoundsImpl_t>) {
                            return insert(std::make_shared<BoundsImpl_t>(args...));
                        }

            private:
                struct BoundComparer{
                    public:
                        bool operator()(const std::shared_ptr<BoundsType_t>& a, 
                                        const std::shared_ptr<BoundsType_t>& b) const {
                            if (a->type() != b->type()) {
                                return static_cast<int>(a->type()) < static_cast<int>(b->type());
                            }
                            const std::vector<double> avalues{a->values()};
                            const std::vector<double> bvalues{b->values()};
                            std::size_t size = avalues.size();
                            for(std::size_t i=0; i<size-1; ++i) {
                                if(std::abs(avalues[i]- bvalues[i]) > std::numeric_limits<double>::epsilon()){
                                    return avalues[i] < bvalues[i];
                                }
                            }
                            return avalues[size-1] < bvalues[size-1];
                        }

            };
            std::set<std::shared_ptr<BoundsType_t>, BoundComparer> m_boundSet{};
        };
    
    using SurfaceBoundFactory = BoundFactory<SurfaceBounds>;
    using VolumeBoundFactory = BoundFactory<VolumeBounds>;
   

}