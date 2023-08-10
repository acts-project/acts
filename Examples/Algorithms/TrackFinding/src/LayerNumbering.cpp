
void ActsExamples::FillInputVector(ActsExamples::SeedingFTFAlgorithm::Config cfg) m_cfg(std::move(cfg)) {
 ///////Hough code ///////////
 ///(sourceLinks = m_cfg.inputSourceLinks)
 ///const auto& measurements = m_inputMeasurements(ctx);
  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) {
    // select volume/layer depending on what is set in the geometry id
   
    auto range = selectLowestNonZeroGeometryObject(sourceLinks, geoId); 

    auto groupedByModule = makeGroupBy(range, detail::GeometryIdGetter()); /// Iterate over groups of elements belonging to each module/ sensitive surface.

    for (auto [moduleGeoId, moduleSourceLinks] : groupedByModule) { 
      // find corresponding surface
      const Acts::Surface* surface =
          m_cfg.trackingGeometry->findSurface(moduleGeoId);
      if (surface == nullptr) {
        ACTS_ERROR("Could not find surface " << moduleGeoId);
        return;
      }

      for (auto& sourceLink : moduleSourceLinks) {

        auto [localPos, localCov] = std::visit(
            [](const auto& meas) {
              auto expander = meas.expander();
              Acts::BoundVector par = expander * meas.parameters();
              Acts::BoundSymMatrix cov =
                  expander * meas.covariance() * expander.transpose();
              // extract local position
              Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
              // extract local position covariance.
              Acts::SymMatrix2 lcov =
                  cov.block<2, 2>(Acts::eBoundLoc0, Acts::eBoundLoc0);
              return std::make_pair(lpar, lcov);
            },
            measurements[sourceLink.index()]);

        // transform local position to global coordinates
        Acts::Vector3 globalFakeMom(1, 1, 1);
        Acts::Vector3 globalPos = surface->localToGlobal(ctx.geoContext, localPos, globalFakeMom); //is this the center? 


        double r = std::hypot(globalPos[Acts::ePos0], globalPos[Acts::ePos1]);
        double phi = std::atan2(globalPos[Acts::ePos1], globalPos[Acts::ePos0]);
        double z = globalPos[Acts::ePos2];
        ResultUnsigned hitlayer = m_cfg.layerIDFinder(r); //is this the layer ACTS layer?? (.value()) //cant find where this is defined 
        if (hitlayer.ok()) {
          std::vector<Index> index;
          index.push_back(sourceLink.index());
          auto meas = std::shared_ptr<HoughMeasurementStruct>(
              new HoughMeasurementStruct(hitlayer.value(), phi, r, z, index, HoughHitType::MEASUREMENT));
          houghMeasurementStructs.push_back(meas);
        }
      }
    }
  }
}
