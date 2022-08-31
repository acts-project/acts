Track Seeding
==============

.. attention::
   This section is largely **outdated** and will be replaced in the future.

To reduce the time needed to reconstruct particle tracks, a track seed
(henceforth: seed) is created which serves as initial direction for the track
reconstruction algorithm (henceforth: the tracking). The tracking then tries to
find all measurements belonging to a single particle in this direction in order
to reconstruct the track. This means, if no seed exists for a particle, this
particle will not be reconstructed. On the other hand, finding too many seeds
which either correspond to a particle for which another seed already exists or
which do not correspond to a particle at all increases the time needed for
tracking.

A good seeding algorithm therefore has the following properties:

#. It finds at least one seed for each particle that should be found
#. It doesn’t find many seeds which do NOT correspond to particles
#. It doesn’t find many seeds per particle

The most typical way to create seeds is to combine measurements. In a homogeneous magnetic field, 3 measurements perfectly describe the helical path of a charged particle. One such triplet of measurements would then constitute a seed and define in close bounds where the tracking needs to look for additional measurements to create a track spanning the whole detector. The difficulty is in choosing the correct measurements, as a helix can be fitted through any 3 measurements in a collision event with potentially tens of thousands of measurements. Therefore, a number of constraints or “cuts” are defined to reduce the number of candidates. Cuts may define where particles originate or the range of energy of particles to be found or otherwise restrict the combination of measurements for seed creation.

Acts Implementation
-------------------

The seeding implementation in Core/include/Acts/Seeding/ is based on the ATLAS track seeding. It was rewritten with a focus on parallelism and maintainability and as detector agnostic as possible, only assuming a (near) homogeneous magnetic field with particles originating from the central detector region. Cuts are configurable and can be plugged in as algorithm which is called by the seeding. The seeding works on measurements or “SpacePoints” (SP), which need to provide x,y,z coordinates with the z axis being along the magnetic field, and x and y. The interaction region must be close to :math:`x=y=0`, such that the interaction region has a smaller detector radius :math:`r = \sqrt{x^2+y^2}` than the measurements closest to the interaction region, see also :numref:`3dim`.

.. figure:: ../figures/seeding/3Dcoordinates.svg
   :name: 3dim
   :align: center
   :width: 600

   Example detector with x,y,z coordinate system. While x and y coordinates must
   be 0 close to the interaction region, z can be arbitrary. The magnetic field
   must be along the z-axis such that charged particles are deflected in \phi
   around the z-axis.

.. figure:: ../figures/seeding/x-yCoordinates.svg
   :name: xy
   :align: center
   :width: 500

   The same detector as above but in x/y coordinate system. Helices show up as
   circles in this projection, as shown by the purple particle path. Three
   measurements define a circle, which is used to calculate radius
   (corresponding to the particle energy), curvature direction in :math:`\phi`
   corresponding to the particle charge, and impact parameters.

.. figure:: ../figures/seeding/r-zCoordinates.svg
   :name: rz
   :align: center
   :width: 500

   The same detector as above but mapped onto an r/z coordinate system. This
   projection is used to calculate the pseudorapidity :math:`\eta` of a seed (in
   the code calculated in :math:`cot \theta` for speed), e.g. to test if two
   combinations (SP bottom, SP middle) and (SP middle, SP top) have similar
   pseudorapidity and are therefore compatible with the same particle track.
   
The search for SPs is performed for each `(\phi, z)` bin of a fully configurable
grid of the detector containing the SPs. The central SP is taken from the
`(\phi, z)` region, then the top and bottom SP are searched in neighboring bins
at the top and bottom layers of the detector, respectively. The number of
neighboring bins used in the SP search for each `(\phi, z)` bin can be defined
in the `z` direction with the vectors `zBinNeighborsTop` and `zBinNeighborsBottom`
separately for bottom and top SPs, and in the `\phi` direction with `numPhiNeighbors`.

Three iterators over SP need to be passed to the public createSeedsForGroup
function in the Seedfinder class. The seedfinder will then attempt to create
seeds, with each seed containing exactly one SP returned by each of the three
iterators. 

- SPs from the first iterator are always used as measurement of a seed with the
  smallest detector radius r, 
- SPs from the second iterator are only used as measurement of a seed with r
  between the r of the first and the third iterator
- SPs from the third iterator are always used as measurement with the largest r
  in a seed.

.. warning::
   Note that the seeding algorithm breaks down for particles with a particle
   track whose helix diameter is smaller than the detector radius until which
   seeds are to be created. This is due to ordering assumptions of SP
   locations as well as due to approximations which become inaccurate for
   lower energy particles.

The createSeedsForGroup function then iterates over middle SP, and within this
loop separately iterates once over bottom SP and once over top SP. Within each
of the nested loops, bottom SP - middle SP respectively middle SP - top SP are
tested for compatibility by applying configurable cuts that can be tested with
two SP only (pseudorapidity, origin along z-axis, distance in r between SP,
compatibility with interaction point).

If both compatible bottom and top SP have been found, test each bottom SP,
middle SP, top SP triplet combination in a triple nested loop. A major part of
this is the calculation of the helix circle. In order to perform calculations
only once, the circle calculation is spread out over the three loops.


.. code-block:: cpp

	for (auto spM : middleSPs) {
   
		... // compatibility cuts between SP duplets

		state.linCircleBottom.clear();
    state.linCircleTop.clear();

		// transform a vector of spacepoints to u-v space circles with respect to a given middle spacepoint
    transformCoordinates(state.compatBottomSP, *spM, true,
                         state.linCircleBottom);
    transformCoordinates(state.compatTopSP, *spM, false, state.linCircleTop);

    state.topSpVec.clear();
    state.curvatures.clear();
    state.impactParameters.clear();
    state.seedsPerSpM.clear();

    size_t numBotSP = state.compatBottomSP.size();
    size_t numTopSP = state.compatTopSP.size();

    int numQualitySeeds = 0;
    int numSeeds = 0;

    size_t t0 = 0;

    for (size_t b = 0; b < numBotSP; b++) {
      auto lb = state.linCircleBottom[b];
      float Zob = lb.Zo;
      float cotThetaB = lb.cotTheta;
      float Vb = lb.V;
      float Ub = lb.U;
      float ErB = lb.Er;
      float iDeltaRB = lb.iDeltaR;

      // 1+(cot^2(theta)) = 1/sin^2(theta)
      float iSinTheta2 = (1. + cotThetaB * cotThetaB);
      // calculate max scattering for min momentum at the seed's theta angle
      // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
      // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
      // scattering
      // but to avoid trig functions we approximate cot by scaling by
      // 1/sin^4(theta)
      // resolving with pT to p scaling --> only divide by sin^2(theta)
      // max approximation error for allowed scattering angles of 0.04 rad at
      // eta=infinity: ~8.5%
      float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
      // multiply the squared sigma onto the squared scattering
      scatteringInRegion2 *=
          m_config.sigmaScattering * m_config.sigmaScattering;

      float sinTheta = 1 / std::sqrt(iSinTheta2);
      float cosTheta = cotThetaB * sinTheta;

      // clear all vectors used in each inner for loop
      state.topSpVec.clear();
      state.curvatures.clear();
      state.impactParameters.clear();
      for (size_t t = t0; t < numTopSP; t++) {
        auto lt = state.linCircleTop[t];

				...	\\ more code
				
				float dU;
        float A;
        float S2;
        float B;
        float B2;

        if (m_config.useDetailedDoubleMeasurementInfo) {
          dU = ut - ub;
          // protects against division by 0
          if (dU == 0.) {
            continue;
          }
          A = (vt - vb) / dU;
          S2 = 1. + A * A;
          B = vb - A * ub;
          B2 = B * B;
        } else {
          dU = lt.U - Ub;
          // protects against division by 0
          if (dU == 0.) {
            continue;
          }
          // A and B are evaluated as a function of the circumference parameters
          // x_0 and y_0
          A = (lt.V - Vb) / dU;
          S2 = 1. + A * A;
          B = Vb - A * Ub;
          B2 = B * B;
        }

        // sqrt(S2)/B = 2 * helixradius
        // calculated radius must not be smaller than minimum radius
        if (S2 < B2 * m_config.minHelixDiameter2) {
          continue;
        }

To calculate the helix circle in the x-y plane, the x,y coordinates are
transformed into a U/V plane in order to calculate the circle with a linear
instead of a quadratic equation for speed. From the helix circle, particle
energy and impact parameters can be estimated.

The scattering calculation is also spread over the nested loops to avoid
redoing calculations. First, the maximum allowed scattering at the configured
minimum transverse momentum (pT) cut is calculated and scaled by the
pseudorapidity of the bottomSP-middleSP duplet to get the minimum momentum of
the duplet. This duplet's pseudorapidity is used for later calculation of the
scattering for the triplet as well.

.. code-block:: cpp

   // 1+(cot^2(theta)) = 1/sin^2(theta)
   float iSinTheta2 = (1. + cotThetaB * cotThetaB);
   // calculate max scattering for min momentum at the seed's theta angle
   // scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
   // accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
   // scattering
   // but to avoid trig functions we approximate cot by scaling by
   // 1/sin^4(theta)
   // resolving with pT to p scaling --> only divide by sin^2(theta)
   // max approximation error for allowed scattering angles of 0.04 rad at
   // eta=infinity: ~8.5%
   float scatteringInRegion2 = m_config.maxScatteringAngle2 * iSinTheta2;
   // multiply the squared sigma onto the squared scattering
   scatteringInRegion2 *=
       m_config.sigmaScattering * m_config.sigmaScattering;
       
The minimum scattering term ('scatteringInRegion2') is calculated from
'sigmaScattering', the configurable number of sigmas of scattering angle
to be considered, and 'maxScatteringAngle2', which is evaluated from the
Lynch & Dahl correction of the Highland equation assuming the lowest
allowed pT. The parameters of the Highland equation are fully configurable.

The following code block calculates if the triplet forms a nearly straight line
in the r/z plane (see :numref:`rz`) as the particle path in the r/z plane is
unaffected by the magnetic field [#f1]_. This is split in two (may be revised
in the future); the first test occurs before the calculation of the helix
circle. Therefore, the deviation from a straight line is compared to the
maximum allowed scattering at minimum pT scaled by the forward angle (as
calculated above). Both the check against min pT and the check against the
calculated pT (discussed further below) take the correlated measurement
uncertainty into account.

.. code-block:: cpp

	// add errors of spB-spM and spM-spT pairs and add the correlation term
	// for errors on spM
	float error2 = lt.Er + ErB +
								 2 * (cotThetaAvg2 * varianceRM + varianceZM) * iDeltaRB *
										 lt.iDeltaR;

	float deltaCotTheta = cotThetaB - cotThetaT;
	float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
	// Apply a cut on the compatibility between the r-z slope of the two
	// seed segments. This is done by comparing the squared difference
	// between slopes, and comparing to the squared uncertainty in this
	// difference - we keep a seed if the difference is compatible within
	// the assumed uncertainties. The uncertainties get contribution from
	// the  space-point-related squared error (error2) and a scattering term
	// calculated assuming the minimum pt we expect to reconstruct
	// (scatteringInRegion2). This assumes gaussian error propagation which
	// allows just adding the two errors if they are uncorrelated (which is
	// fair for scattering and measurement uncertainties)
	if (deltaCotTheta2 > (error2 + scatteringInRegion2)) {
		// additional cut to skip top SPs when producing triplets
		if (m_config.skipPreviousTopSP) {
			// break if cotTheta from bottom SP < cotTheta from top SP because
			// the SP are sorted by cotTheta
			if (cotThetaB - cotThetaT < 0) {
				break;
			}
			t0 = t + 1;
		}
		continue;
	}

Now the check for scattering with calculated particle momentum. Momentum is
calculated from the pT and the pseudorapidity. This must be :math:`\geq` the
lower pT cut, and therefore allows :math:`\leq` scattering compared to the
previous check, as the scattering decreases linearly with particle energy

.. code-block:: cpp

	// refinement of the cut on the compatibility between the r-z slope of
	// the two seed segments using a scattering term scaled by the actual
	// measured pT (p2scatterSigma)
	float iHelixDiameter2 = B2 / S2;
	// calculate scattering for p(T) calculated from seed curvature
	float pT2scatterSigma = iHelixDiameter2 * m_config.sigmapT2perRadius;
	// if pT > maxPtScattering, calculate allowed scattering angle using
	// maxPtScattering instead of pt.
	float pT = m_config.pTPerHelixRadius * std::sqrt(S2 / B2) / 2.;
	if (pT > m_config.maxPtScattering) {
		float pTscatterSigma =
				(m_config.highland / m_config.maxPtScattering) *
				m_config.sigmaScattering;
		pT2scatterSigma = pTscatterSigma * pTscatterSigma;
	}
	// convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
	// from rad to deltaCotTheta
	float p2scatterSigma = pT2scatterSigma * iSinTheta2;
	// if deltaTheta larger than allowed scattering for calculated pT, skip
	if (deltaCotTheta2 > (error2 + p2scatterSigma)) {
		if (m_config.skipPreviousTopSP) {
			if (cotThetaB - cotThetaT < 0) {
				break;
			}
			t0 = t;
		}
		continue;
	}

The last cut applied in this function is on the so-called impact parameters,
which is the distance of the perigee of a track from the interaction region in
mm of detector radius. It is calculated and cut on before storing all top SP
compatible with both the current middle SP and current bottom SP.

.. code-block:: cpp

	// A and B allow calculation of impact params in U/V plane with linear
	// function
	// (in contrast to having to solve a quadratic function in x/y plane)
	float Im;
	if (m_config.useDetailedDoubleMeasurementInfo == false) {
		Im = std::abs((A - B * rM) * rM);
	} else {
		Im = std::abs((A - B * rMxy) * rMxy);
	}

	if (Im <= m_config.impactMax) {
		state.topSpVec.push_back(state.compatTopSP[t]);
		// inverse diameter is signed depending if the curvature is
		// positive/negative in phi
		state.curvatures.push_back(B / std::sqrt(S2));
		state.impactParameters.push_back(Im);

		// evaluate eta and pT of the seed
		float cotThetaAvg = std::sqrt(cotThetaAvg2);
		float theta = std::atan(1. / cotThetaAvg);
		float eta = -std::log(std::tan(0.5 * theta));
		state.etaVec.push_back(eta);
		state.ptVec.push_back(pT);
	}

The bottom SP and middle SP as well as the collection of top SP is passed to
SeedFilter::filterSeeds_2SpFixed, whose collected output for the current middle
SP with all compatible bottom SP and top SP is then passed to
SeedFilter::filterSeeds_1SpFixed.

SeedFilter::filterSeeds_2SpFixed
--------------------------------

This function assigns a weight (which should correspond to the likelihood that
a seed is good) to all seeds and calls the detector specific cuts to apply a
hard cut or modify the weight. The weight is a “soft cut”, in that it is only
used to discard tracks if many seeds are created for the same middle SP in
SeedFilter::filterSeeds_1SpFixed. This process is important to improving computational
performance and the quality of the final track collections by rejecting lower-quality seeds.

The weight can be influenced by:

#. The transverse (`d_{0}`) and longitudinal (`z_{0}`) impact parameters (the higher the distance the worse)
#. The number of seeds which may belong to the same particle track (`N_{t}`)
#. Optional detector specific cuts.

The transverse impact parameter is multiplied by the configured factor and subtracted from
the weight, as seeds with higher impact parameters are assumed to be less
likely to stem from a particle than another seed using the same middle SP with
smaller impact parameters. The longitudinal impact parameter is subtracted from
the weight if configured.

The number of seeds only differing in top SP which have similar helix radius
and the same sign (i.e. the same charge) is used to increase the weight, as it
means that more than three measurements that may be from the same particle have
been found. The measurements must have a minimum distance in detector radius,
such that measurements from the same layer cannot be counted towards the
increased weight. The number of found compatible seeds is multiplied by a
configured factor and added to the weight.

The optional detector specific cuts can use the weight from 1. and 2. and the
three SP to apply a hard cut or change the weight of a seed.

This function also includes a fully configurable seed confirmation step that, when enabled
('seedConfirmation=True'), classifies higher quality seeds as "quality confined" seeds if
they fall within a predefined range of parameters (`d_{0}`, `z_{0}` and `N_{t}`) that also
depends on the region of the detector (i.e., forward or central region). If the seed is not
classified as "quality confined" seed, it will only be accepted if its weight is greater
than a certain threshold and no other high quality seed has been found.

.. code-block:: cpp

    int deltaSeedConf;
    if (m_cfg.seedConfirmation) {
      // seed confirmation cuts - keep seeds if they have specific values of
      // impact parameter, z-origin and number of compatible seeds inside a
      // pre-defined range that also depends on the region of the detector (i.e.
      // forward or central region) defined by SeedConfirmationRange
      deltaSeedConf = compatibleSeedR.size() + 1 - nTopSeedConf;
      if (deltaSeedConf < 0 || (numQualitySeeds and deltaSeedConf == 0)) {
        continue;
      }
      bool seedRangeCuts = bottomSP.radius() < m_cfg.seedConfMinBottomRadius ||
                           std::abs(zOrigin) > m_cfg.seedConfMaxZOrigin;
      if (seedRangeCuts and deltaSeedConf == 0 and
          impact > m_cfg.minImpactSeedConf) {
        continue;
      }

The seed confirmation also sets a limit on the number of seeds produced for each middle SP,
which retains only the higher quality seeds. If this limit is exceeded, the algorithm
checks if there is any low-quality seed in the seed container of this middle SP that can be removed.

SeedFilter::filterSeeds_1SpFixed
--------------------------------

This function allows the detector specific cuts to filter on the basis of all
seeds with a common middle SP and limits the number of seeds per middle SP to
the configured limit. It sorts the seeds by weight and, to achieve a
well-defined ordering in the rare case weights are equal, sorts them by
location. The ordering by location is only done to make sure reimplementations
(such as the GPU code) are comparable and return the bitwise exactly same
result.

When a seed is accepted and seed confirmation is enabled, the weight of that seed
is assigned to each of its SPs. Each SP will hold the weight of the best seed that
includes that SP. This information is used in the selection of the next seeds:
The seed is kept only if its weight is greater than the weight of at least one of
its SP components.


Footnotes
---------

.. [#f1] approximately, this is one of the reasons the algorithm breaks down for low energy particles.

