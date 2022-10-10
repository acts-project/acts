void Acts::Experimental::connectCylindricalVolumes(
    const GeometryContext& gctx,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const CylindricalDetectorHelperOptions& options) {
  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("DetectorVolumeConverter", options.logLevel));

  // Consistency checking block -------------
  //
  // First check:
  // - type of bounds
  // - alginment (z-axis aligned)
  //
  Vector3 refZaxis(0., 0., 1);
  std::array<Vector3, 2> refZspan;

  // Pairwise check
  Vector3 cogPW(0., 0., 0.);
  ActsScalar accPhi = 0.;

  std::array<ActsScalar, 2> refRadii;
  std::array<ActsScalar, 2> refPhi = {M_PI, 0.};

  std::unordered_set<BinningValue> attachmentOptions = {binR, binZ, binPhi};

  for (const auto [i, v] : enumerate(volumes)) {
    const auto& vBounds = v->volumeBounds();
    const auto center = v->transform(gctx).translation();
    cogPW += center;
    // Type checking on the bounds
    if (vBounds.type() != Acts::VolumeBounds::eCylinder) {
      throw std::invalid_argument(
          "CylindricalVolumeHelper: volumes are not of cylindrical shape.");
    }
    accPhi +=
        2 * vBounds.values()[CylinderVolumeBounds::BoundValues::eHalfPhiSector];

    ActsScalar minR =
        vBounds.values()[CylinderVolumeBounds::BoundValues::eMinR];
    ActsScalar maxR =
        vBounds.values()[CylinderVolumeBounds::BoundValues::eMaxR];

    // Axis check
    if (i == 0u) {
      refRadii = {minR, maxR};
      refZaxis = v->transform(gctx).rotation().matrix().col(2);
      refZspan[0u] =
          center -
          vBounds.values()[CylinderVolumeBounds::BoundValues::eHalfLengthZ] *
              refZaxis;
      refZspan[1u] =
          center +
          vBounds.values()[CylinderVolumeBounds::BoundValues::eHalfLengthZ] *
              refZaxis;

    } else {
      if (not refZaxis.isApprox(
              v->transform(gctx).rotation().matrix().col(2))) {
        throw std::invalid_argument(
            "CylindricalVolumeHelper: cylindrical volumes are not "
            "aligned.");
      }
      refZspan[1u] =
          center +
          vBounds.values()[CylinderVolumeBounds::BoundValues::eHalfLengthZ] *
              refZaxis;

      // Update new center of gravity
      cogPW /= 2u;
      // Check if the new center is consistent, otherwise erase the binZ option
      if (cogPW.isApprox(center)) {
        if (attachmentOptions.count(binZ) > 0u) {
          attachmentOptions.erase(binZ);
          ACTS_VERBOSE("Binning in Z option excluded.")
        }
      } else if (attachmentOptions.count(binR) > 0u) {
        attachmentOptions.erase(binR);
        ACTS_VERBOSE("Binning in R option excluded.")
      }
      if (maxR <= refRadii[1] and minR >= refRadii[0] and
          attachmentOptions.count(binR) > 0u) {
        attachmentOptions.erase(binR);
      }
      // Check accumulated opening angle, otherwise erase the binPhi option
      if (accPhi > 2 * M_PI + 5 * std::numeric_limits<ActsScalar>::epsilon() and
          attachmentOptions.count(binPhi) > 0u) {
        attachmentOptions.erase(binPhi);
        ACTS_VERBOSE("Binning in phi option excluded.")
      }
    }
  }

  if (attachmentOptions.size() != 1u) {
    throw std::invalid_argument(
        "CylinderVolumeHelper: attachment options could not be determined.");
  }

  BinningValue aOption = *(attachmentOptions.begin());

  std::vector<std::shared_ptr<DetectorVolume>> internalVolumes = volumes;

  // Consistency checking under attachment options known, the volumes are
  // here assumed to be ordered in attachment direction

  // The pValues are the parameters in the binning direction
  std::vector<std::array<ActsScalar, 2>> pValues = {};
  // Referencing
  Vector3 lastCenter(0., 0., 0.);
  const VolumeBounds* lastBounds =
      &(internalVolumes.begin()->get()->volumeBounds());
  switch (aOption) {
    // Halflength in z needs to be the same, phi opening compatible
    case binR: {
      for (const auto [i, v] : enumerate(internalVolumes)) {
        if (i > 0) {
          if (std::abs(lastBounds->values()
                           [CylinderVolumeBounds::BoundValues::eHalfLengthZ] -
                       v->volumeBounds().values()
                           [CylinderVolumeBounds::BoundValues::eHalfLengthZ]) >
              options.glueTolerance) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: r-binned volumes need to have the same "
                "half lenghts.");
          }
          if (std::abs(
                  lastBounds->values()
                      [CylinderVolumeBounds::BoundValues::eHalfPhiSector] -
                  v->volumeBounds().values()
                      [CylinderVolumeBounds::BoundValues::eHalfPhiSector]) >
              options.glueTolerance) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: r-binned volumes need to have the same "
                "phi openings.");
          }
        }
        // The principle values
        pValues.push_back(
            {v->volumeBounds()
                 .values()[CylinderVolumeBounds::BoundValues::eMinR],
             v->volumeBounds()
                 .values()[CylinderVolumeBounds::BoundValues::eMaxR]});
        lastBounds = &(v->volumeBounds());
      }
    } break;

    // Radii need to be the same and, phi opening angles compatible
    case binZ: {
      for (const auto [i, v] : enumerate(internalVolumes)) {
        if (i == 0) {
          pValues.push_back(
              {-(v->volumeBounds().values()
                     [CylinderVolumeBounds::BoundValues::eHalfLengthZ]),
               v->volumeBounds()
                   .values()[CylinderVolumeBounds::BoundValues::eHalfLengthZ]});

        } else {
          // Distance between first and last center, reference is first center
          ActsScalar dist = (v->transform(gctx).translation() -
                             internalVolumes[0u]->transform(gctx).translation())
                                .norm();
          pValues.push_back(
              {dist - (v->volumeBounds().values()
                           [CylinderVolumeBounds::BoundValues::eHalfLengthZ]),
               dist + v->volumeBounds().values()
                          [CylinderVolumeBounds::BoundValues::eHalfLengthZ]});
          // R values need to be consistent
          if ((std::abs(
                   lastBounds
                       ->values()[CylinderVolumeBounds::BoundValues::eMinR] -
                   v->volumeBounds()
                       .values()[CylinderVolumeBounds::BoundValues::eMinR]) >
               options.glueTolerance) or
              (std::abs(
                   lastBounds
                       ->values()[CylinderVolumeBounds::BoundValues::eMaxR] -
                   v->volumeBounds()
                       .values()[CylinderVolumeBounds::BoundValues::eMaxR]) >
               options.glueTolerance)) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: z-binned volumes need to have the "
                "same r values.");
          }
          // Phi sectors needs to be consistent
          if (std::abs(
                  lastBounds->values()
                      [CylinderVolumeBounds::BoundValues::eHalfPhiSector] -
                  v->volumeBounds().values()
                      [CylinderVolumeBounds::BoundValues::eHalfPhiSector]) >
              options.glueTolerance) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: z-binned volumes need to have the same "
                "phi openings.");
          }
        }
        lastBounds = &(v->volumeBounds());
        lastCenter = v->transform(gctx).translation();
      }
    } break;
    case binPhi: {
      for (const auto [i, v] : enumerate(internalVolumes)) {
        if (i > 0) {
          if ((std::abs(lastBounds->values()
                            [CylinderVolumeBounds::BoundValues::eHalfLengthZ] -
                        v->volumeBounds().values()
                            [CylinderVolumeBounds::BoundValues::eHalfLengthZ]) >
               options.glueTolerance)) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: phi-binned volumes need to have the "
                "same half lengths.");
          }
          if ((std::abs(
                   lastBounds
                       ->values()[CylinderVolumeBounds::BoundValues::eMinR] -
                   v->volumeBounds()
                       .values()[CylinderVolumeBounds::BoundValues::eMinR]) >
               options.glueTolerance) or
              (std::abs(
                   lastBounds
                       ->values()[CylinderVolumeBounds::BoundValues::eMaxR] -
                   v->volumeBounds()
                       .values()[CylinderVolumeBounds::BoundValues::eMaxR]) >
               options.glueTolerance)) {
            throw std::invalid_argument(
                "CylinderVolumeHelper: phi-binned volumes need to have the "
                "same r values.");
          }
          lastBounds = &(v->volumeBounds());
          lastCenter = v->transform(gctx).translation();
        }
      }
      break;
      default:
        break;
    }
  }

  auto actualHandling = CylindricalDetectorHelperOptions::Handling::eStrict;

  // Consistency checking of p values
  std::vector<ActsScalar> pBoundaries = {};  // these are the binning boundaries
  std::vector<ActsScalar> pReferences = {};  // these are the binning references
  std::array<ActsScalar, 2> lValues;
  std::array<ActsScalar, 2> mValues;
  bool resizeOrFill = false;
  for (auto [ip, p] : enumerate(pValues)) {
    if (ip > 0u) {
      if (lValues[1] > p[0] + std::numeric_limits<ActsScalar>::epsilon()) {
        throw std::invalid_argument(
            "CylinderVolumeHelper: overlapping values!");
      }
      if (std::abs(lValues[1] - p[0]) >
          std::numeric_limits<ActsScalar>::epsilon()) {
        ACTS_DEBUG("Resizing of Filling necessary as " << lValues[1]
                                                       << " =/= " << p[0]);
        resizeOrFill = true;
        // Resizing: take the middle value
        if (options.handling ==
            CylindricalDetectorHelperOptions::Handling::eResize) {
          pBoundaries.push_back(0.5 * (lValues[1] + p[0]));
        } else if (options.handling ==
                   CylindricalDetectorHelperOptions::Handling::eFill) {
          // Fill: make a gap volume
          pBoundaries.push_back(lValues[1]);
          pBoundaries.push_back(p[0]);
        }
      } else {
        pBoundaries.push_back(p[0]);
      }
    }
    lValues = p;
    mValues = {std::min(mValues[0], p[0]), std::max(mValues[1], p[1])};
  }
  // All boundaries filled
  pBoundaries.push_back(mValues[0]);
  pBoundaries.push_back(mValues[1]);
  std::sort(pBoundaries.begin(), pBoundaries.end());

  // Bail out if exact handling is espected, but not possible
  if (not resizeOrFill) {
    ACTS_DEBUG("Exact attachment detected, no resizing necessary");
    actualHandling = CylindricalDetectorHelperOptions::Handling::eStrict;
  } else if (options.handling ==
             CylindricalDetectorHelperOptions::Handling::eStrict) {
    throw std::invalid_argument(
        "CylinderVolumeHelper: resizing necessary, but strict handling "
        "chosen.");
  }

  ActsScalar pReference = 0.;
  // If handling is not yet strict, resize or insert until strict
  if (actualHandling != CylindricalDetectorHelperOptions::Handling::eStrict) {
    /// Hada hada hada
  }

  switch (aOption) {
    case binR: {
      ACTS_VERBOSE("Attach with r binning option.");
      connectVolumesInR(gctx, internalVolumes, pBoundaries);
    } break;
    case binZ: {
      ACTS_VERBOSE("Attach with z binning option.");
      connectVolumesInZ(gctx, internalVolumes, pBoundaries);
    } break;
    case binPhi: {
      ACTS_VERBOSE("Attach with phi binning option.");
      connectVolumesInPhi(gctx, internalVolumes, pBoundaries);
    } break;

    default:
      break;
  }

  // Consistency checking end
}