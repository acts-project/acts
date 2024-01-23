// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/container/static_vector.hpp>

#include "Acts/Detector/Detector.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Detector/DetectorVolumeBuilder.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdaters.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/EventData/detail/TestSourceLink.hpp"
#include "Acts/EventData/detail/GenerateParameters.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderOptions.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/SurfaceContainer.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"

#include "ActsExamples/Utilities/GroupBy.hpp"

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <functional>

namespace Acts {
namespace Test {
    // make detector
    using namespace UnitLiterals;
    using namespace detail::Test;
    using namespace ActsExamples;
    using BuilderFunction = std::function<ActsExamples::SimSpacePoint(
            Acts::Vector3, std::optional<double>, Acts::Vector2,
            std::optional<double>, boost::container::static_vector<SourceLink, 2>)>;
    // make tracking geometry
    GeometryContext gctx = GeometryContext();
//        BOOST_AUTO_TEST_SUITE(SurfaceContainer)

    BOOST_AUTO_TEST_CASE(SurfaceContainerTest){
            // Make basic Objects
            unsigned numSurfaces = 6;
            auto rBounds = std::make_shared<const RectangleBounds>(3,4);
            std::vector<std::shared_ptr<Acts::PlaneSurface>> surfaces;
            std::vector<Transform3> transformations;
            std::array<LayerPtr, 6> layers{};
            LayerArrayCreator::Config lacConfig;
            LayerArrayCreator layArrCreator(
            lacConfig, getDefaultLogger("LayerArrayCreator", Logging::INFO));
            Vector2 Zero2{ 0., 0. };
            Vector3 Zero3{ 0., 0., 0. };

            for (unsigned i=0; i<numSurfaces; i++) {
                // Make Surfaces
                long double x = (long double) i;
                Translation3 translation{(x + 1.) * 10._mm, 0.25_m, 0.25_m};
                auto pTransform = Transform3(translation);
                transformations.push_back(pTransform);
                auto planeSurface =
                        Surface::makeShared<PlaneSurface>(pTransform, rBounds);
                surfaces.push_back(planeSurface);
                    BOOST_CHECK_EQUAL(
                            surfaces.at(i)->localToGlobal(gctx, Zero2, Zero3)(0), (x + 1) * 10._mm
                    );
                    BOOST_CHECK_EQUAL(
                            surfaces.at(i)->localToGlobal(gctx, Zero2, Zero3)(1), 0.25_m
                    );
                    BOOST_CHECK_EQUAL(
                            surfaces.at(i)->localToGlobal(gctx, Zero2, Zero3)(2), 0.25_m
                    );

                // Make Layers
                Transform3 trafo(translation * RotationMatrix3::Identity());

                std::unique_ptr<SurfaceArray> surArray(new SurfaceArray(surfaces[i]));

                layers[i] = PlaneLayer::create(trafo, rBounds, std::move(surArray), 0._mm);

                auto mutableSurface = const_cast<PlaneSurface *>(surfaces[i].get());
                mutableSurface->associateLayer(*layers[i]);
                    BOOST_CHECK_EQUAL(
                            ((layers[i]->surfaceRepresentation()).localToGlobal(gctx, Zero2, Zero3)).isApprox(
                                    surfaces.at(i)->localToGlobal(gctx, Zero2, Zero3)), true
                    );
                    BOOST_CHECK_EQUAL(
                            ((layers[i]->surfaceRepresentation()).localToGlobal(gctx, Zero2, Zero3))(0),
                            (x + 1) * 10._mm
                    );
            }

            Transform3 trafoVol(Transform3::Identity());
            trafoVol.translation() = Vector3(0., 0., 0.);

            LayerVector layVec;
            for (unsigned int i = 0; i < 6; i++) {
                layVec.push_back(layers[i]);
            }
            std::unique_ptr<const LayerArray> layArr(
            layArrCreator.layerArray(gctx, layVec, 0., 2._m,
            BinningType::arbitrary, BinningValue::binX));
            auto boundsVolTG =
            std::make_shared<const CuboidVolumeBounds>(1.5_m, 0.5_m, 0.5_m);

            auto trackVolume =
            TrackingVolume::create(trafoVol, boundsVolTG, nullptr,
            std::move(layArr), nullptr, {}, "Volume");

            Transform3 trafoWorld(Transform3::Identity());
            trafoWorld.translation() = Vector3(0., 0., 0.);

            auto worldVol =
            std::make_shared<const CuboidVolumeBounds>(3._m, 0.5_m, 0.5_m);

            std::vector<std::pair<TrackingVolumePtr, Vector3>> tapVec;

            tapVec.push_back(std::make_pair(trackVolume, Zero3));

            std::vector<float> binBoundaries = { -3._m, 0., 3._m };

            BinningData binData(BinningOption::open, BinningValue::binX, binBoundaries);
            std::unique_ptr<const BinUtility> bu(new BinUtility(binData));

            std::shared_ptr<const TrackingVolumeArray> trVolArr(
            new BinnedArrayXD<TrackingVolumePtr>(tapVec, std::move(bu)));

            MutableTrackingVolumePtr mtvpWorld(
            TrackingVolume::create(trafoWorld, worldVol, trVolArr, "World"));
            // Build tracking geometry
            std::shared_ptr<TrackingGeometry> tGeometry(
            new Acts::TrackingGeometry(mtvpWorld));
            // make surfacePtrs from TG
            auto surfacePtrsTG = SurfaceContainer(tGeometry).surfacePtrs();
            BOOST_CHECK_NE(tGeometry, nullptr);
            BOOST_TEST_MESSAGE("PASSED tracking geometry tests" );
            auto portalGenerator = Acts::Experimental::defaultPortalGenerator();
            Acts::Transform3 nominal = Acts::Transform3::Identity();
            BOOST_CHECK_NE(tGeometry, nullptr);

            std::vector<std::shared_ptr<Acts::Surface>> surfacePtrs_shared;
            GeometryIdentifier surfaceId;
            for (unsigned i=0 ; i<numSurfaces ; i++) {
//                auto detElement = std::make_unique<Acts::Test::DetectorElementStub>(
//                        transformations[i],
//                        rBounds, 0.01);
//                auto surface = detElement->surface().getSharedPtr();
                auto detPlaneSurface = Surface::makeShared<PlaneSurface>(transformations[i],rBounds);
                surfaceId.setVolume(2u).setLayer(2u*(i+1u)).setSensitive(1u);
                detPlaneSurface->assignGeometryId(surfaceId);
                surfacePtrs_shared.push_back(detPlaneSurface);
            }

            auto boundsVolDet =
            std::make_shared<Acts::CuboidVolumeBounds>(1.5_m, 0.5_m, 0.5_m);
            //const std::vector<std::shared_ptr<RegularSurface>> surfacePtrs_shared_const(surfacePtrs_shared);

            auto vol = Acts::Experimental::DetectorVolumeFactory::construct(
            portalGenerator, gctx, "VolumeWithSurfaces", nominal,
            std::move(boundsVolDet), surfacePtrs_shared, {},
            Acts::Experimental::tryNoVolumes(),
            Acts::Experimental::tryAllPortalsAndSurfaces());

            Acts::GeometryIdentifier volId;
            volId.setVolume(2u);
            vol->assignGeometryId(volId);

            auto detector = Acts::Experimental::Detector::makeShared(
            "detector", {vol}, Acts::Experimental::tryRootVolumes());

            SurfaceContainer detSC(detector);
            auto surfacePtrsD = detSC.surfacePtrs();
            for (unsigned i=0 ; i<numSurfaces ; i++) {
                auto detTransform = surfacePtrsD.at(i)->localToGlobal(gctx,Zero2,Zero3);
                auto TGTransform =  surfacePtrsTG.at(i)->localToGlobal(gctx,Zero2,Zero3);
                BOOST_CHECK_EQUAL(detTransform.isApprox(TGTransform),true);
            }
    }
    };// namespace Test
} // namcespace Acts
//
//ActsExamples::IndexSourceLinkContainer sourceLinks;
///// geometry and indexing events
//std::vector<::Acts::BoundVariantMeasurement> measOriginal;
//std::vector<double> xValues{ 0.1, 0.2, 0.1, 0.2, 0.1, 0.3 };
//std::vector<double> yValues{ 1.1, 0.2, 1.1, 1.2, 0.1, 1.3 };
////std::vector<unsigned long>
//std::vector<Acts::GeometryIdentifier> gsel;
//for (auto entry = 0; entry < 5; entry++) {
//// Process the correct event
//Acts::GeometryIdentifier geoId;
//geoId.setVolume(2ul).setLayer(2ul * entry).setSensitive(1ul);
//gsel.push_back(geoId);
//
//const Acts::ActsVector<2> par = {
//        xValues[entry] * Acts::UnitConstants::mm,
//        yValues[entry] * Acts::UnitConstants::mm,
//};
//
//Acts::ActsSquareMatrix<2> cov = Acts::ActsSquareMatrix<2>::Identity();
//
//ActsExamples::IndexSourceLink sl(geoId, static_cast<ActsExamples::Index>(entry));
//sourceLinks.insert(sl);
//auto m = Acts::makeMeasurement(Acts::SourceLink{sl}, par, cov, Acts::eBoundLoc0, Acts::eBoundLoc1);
//
//measOriginal.push_back(m);
//
//}
//            auto spBuilderConfig_TG = Acts::SpacePointBuilderConfig();
//            auto spBuilderConfig_SP = Acts::SpacePointBuilderConfig();
//            spBuilderConfig_SP.surfacePtrs = surfacePtrsTG;
//            spBuilderConfig_TG.trackingGeometry = tGeometry;
//
//            std::optional<IndexSourceLink::SurfaceAccessor> m_slSurfaceAccessor;
//            m_slSurfaceAccessor.emplace(
//            IndexSourceLink::SurfaceAccessor{ surfacePtrsTG });
//            spBuilderConfig_TG.slSurfaceAccessor
//            .connect<&IndexSourceLink::SurfaceAccessor::operator()>(
//            &m_slSurfaceAccessor.value());
//
//            BuilderFunction spConstructor =
//            [](Acts::Vector3 pos, std::optional<double>, Acts::Vector2 cov, std::optional<double>,
//            boost::container::static_vector<Acts::SourceLink, 2> slinks)
//            -> ActsExamples::SimSpacePoint {
//                return SimSpacePoint(pos, 0, cov[0], cov[1], 1, std::move(slinks));
//            };
//            Acts::Logging::Level lvl = Acts::Logging::DEBUG;
//
//            Acts::SpacePointBuilder spacePointBuilder_TG = Acts::SpacePointBuilder<SimSpacePoint>(
//            spBuilderConfig_TG, spConstructor,
//            Acts::getDefaultLogger("TestSpacePointBuilder", lvl));
//
//            Acts::SpacePointBuilderOptions spOpt;
//            spOpt.paramCovAccessor =[](Acts::SourceLink slink) { //&measOriginal
//                const auto islink = slink.get<IndexSourceLink>();
//                std::vector<double> xValues{0.1, 0.2, 0.1, 0.2, 0.1, 0.3};
//                std::vector<double> yValues{1.1, 0.2, 1.1, 1.2, 0.1, 1.3};
//                BoundVector par = {
//                        xValues[islink.index()] * Acts::UnitConstants::mm,
//                        yValues[islink.index()] * Acts::UnitConstants::mm,
//                        0.,
//                        0.,
//                        0.,
//                        0.
//                };
//                Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
//                return std::make_pair(par, cov);
//            };
//
//
//            // make points
//            SimSpacePointContainer spacePoints;
//            for (Acts::GeometryIdentifier geoId : gsel) {
//                // select volume/layer depending on what is set in the geometry id
//                if (geoId.sensitive() == 0)
//                    continue;
//                auto range = selectLowestNonZeroGeometryObject(sourceLinks, geoId);
//                //BOOST_CHECK_NE(range.size(),0);
//                // groupByModule only works with geometry containers, not with an
//                // arbitrary range. do the equivalent grouping manually
//                auto groupedByModule = makeGroupBy(range, ActsExamples::detail::GeometryIdGetter());
//                for (auto [moduleGeoId, moduleSourceLinks]: groupedByModule) {
//                    for (auto &sourceLink: moduleSourceLinks) {
//                        spacePointBuilder_TG.buildSpacePoint(
//                                gctx, {Acts::SourceLink{sourceLink}}, spOpt,
//                                std::back_inserter(spacePoints));
//                    }
//                }
//                    BOOST_CHECK_EQUAL(spacePoints.size(),6);
//            }
//            // make spacepoints

