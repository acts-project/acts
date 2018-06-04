#ifndef SeedmakerConfig
#define SeedmakerConfig


namespace Acts{
      // forward declaration to avoid cyclic dependence
      class IBinFinder;
      class ISeedFilter;
      class ICovarianceTool;
      struct Config {
        // BinFinder must return std::vector<Acts::Seeding::Bin> with content of 
        // each bin sorted in r (ascending)
        // TODO: move helper tools to own struct
        std::shared_ptr<IBinFinder> bottomBinFinder;
        std::shared_ptr<IBinFinder> topBinFinder;
        std::shared_ptr<ISeedFilter> seedFilter;
        std::shared_ptr<ICovarianceTool> covarianceTool;

// Seed Cuts
        // lower cutoff for seeds in MeV
        float minPt = 400.;
        float cotThetaMax;
        float deltaRMin;
        float deltaRMax;

        float minPhi = -M_PI;
        float maxPhi = M_PI;

        // the delta for inverse helix radius up to which compared seeds
        // are considered to have a compatible radius. delta of inverse radius
        // leads to this value being the cutoff. unit is 1/mm. default value
        // of 0.00003 leads to all helices with radius>33m to be considered compatible

        // impact parameter in mm
        float impactMax = 20.;

        // how many sigmas of scattering angle should be considered?
        float sigmaScattering = 5;

        // for how many seeds can one SpacePoint be the middle SpacePoint?
        int maxSeedsPerSpM;

// Geometry Settings
        // Detector ROI
        float collisionRegionMin;
        float collisionRegionMax;
        float phiMin = -M_PI;
        float phiMax = M_PI;
        float zMin;
        float zMax;
        // rmax only influences phi bin size
        float rMax;

        // Unit in kiloTesla
        float bFieldInZ = 0.00208;
        // location of beam in x,y plane. 
        // used to moveints
        std::array<float,2> beamPos{{0,0}};

        // average amount of material on the length of a seed. used for scattering.
        // default is 5%
        // TODO: necessary to make amount of material dependent on detector region?
        float radLengthPerSeed = 0.05;
        // alignment uncertainties, used for uncertainties in the non-measurement-plane of the modules 
        // which otherwise would be 0
        // will be added to spacepoint measurement uncertainties (and therefore also multiplied by sigmaError)
        float zAlign = 0;
        float rAlign = 0;
        // used for measurement (+alignment) uncertainties.
        float sigmaError = 5;
      };
} //namespace Acts

#endif //SeedmakerConfig
