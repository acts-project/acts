.. _seeding_core:

Track Seeding
==============

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
   
   Sketch of the detector with 3 layers. The interaction region is supposed to be located along the z axis and have size significantly smaller than the radius of the innermost detector layer.

.. figure:: ../figures/seeding/x-yCoordinates.svg
   :name: xy
   :align: center
   :width: 500
   
   The x-y projection of the detector with the charged particle helical track originating from the centre of the detector. Signals left by passage of the track through the detector layers are marked with green crosses.

.. figure:: ../figures/seeding/r-zCoordinates.svg
   :name: rz
   :align: center
   :width: 500
   
   The r-z projection of the detector with the same charged particle track. The track is depicted with the same colours as on previous figure.
   
The SPs in each detector layer are projected on a rectangular grid of configurable
granularity. The search for seed starts from selecting SP in the middle detector 
layer. Then matching SPs are searched in the inner an outer layers. Grouping of 
the UPs in the aforementioned grid allows to limit the search to neighbouring grid 
cells thus improving significantly algorithm performance. The number of neighboring 
bins used in the SP search for each `(\phi, z)` bin can be defined
in the `z` direction with the vectors `zBinNeighborsTop` and `zBinNeighborsBottom`
separately for bottom and top SPs, and in the `\phi` direction with `numPhiNeighbors`.

The method to create the seed is `createSeedsForGroup`. It receives three iterators 
over SPs constructed from detector layers of increasing radii. The seedfinder will 
then attempt to create seeds, with each seed containing exactly one SP returned by 
each of the three iterators. 

.. warning::
   Note that the seeding algorithm breaks down for particles with a particle
   track whose helix diameter is smaller than the detector radius until which
   seeds are to be created. This is due to ordering assumptions of SP
   locations as well as due to approximations which become inaccurate for
   lower energy particles.

The `createSeedsForGroup` function then iterates over SPs in the middle layer
(2nd iterator), and within this loop separately iterates once over bottom SP 
and once over top SP. Within each of the nested loops, SP pairs are tested for
compatibility by applying a set of configurable cuts that can be tested with
two SP only (pseudorapidity, origin along z-axis, distance in r between SP,
compatibility with interaction point).

If both compatible bottom and top SP have been found, test each bottom SP,
middle SP, top SP triplet combination in a triple nested loop. A major part of
this is the calculation of the helix circle. In order to perform calculations
only once, the circle calculation is spread out over the three loops.


.. code-block:: cpp

	for (auto spM : middleSPs) {
   
		// ... // compatibility cuts between SP duplets

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

				// ...	\\ more code
				
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

From the helix circle, particle energy and impact parameters can be estimated.
To calculate the helix circle in the :math:`x/y` plane, the x,y coordinates are
transformed into a :math:`u/v` plane in order to calculate the circle with a linear equation
instead of a quadratic equation for speed. The conformal transformation is given by:

$$
u = \\frac{x}{x^2+y^2}, \\quad \\quad v = \\frac{y}{x^2+y^2}
$$

Where the circle containing the three SPs are transformed into a line with equation :math:`v = Au + B`


The scattering calculation is also spread over the nested loops to avoid
redoing calculations. First, the maximum allowed scattering at the configured
minimum transverse momentum (:math:`p_{T}`) cut is calculated and scaled by the
pseudorapidity of the duplet formed by one SP from bottom layer and one SP from middle layer to get the minimum momentum of
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
       
The minimum scattering term (`scatteringInRegion2`) is calculated from
`sigmaScattering`, the configurable number of sigmas of scattering angle
to be considered, and `maxScatteringAngle2`, which is evaluated from the
Lynch & Dahl correction of the Highland equation assuming the lowest
allowed :math:`p_{T}`. The parameters of the Highland equation are fully configurable.

The following code block checks if the triplet forms a nearly straight line
in the :math:`r/z` plane (see :numref:`rz`) as the particle path in the :math:`r/z` plane is
unaffected by the magnetic field [#f1]_. This is split in two parts (may be revised
in the future); the first test occurs before the calculation of the helix
circle. Therefore, the deviation from a straight line is compared to the
maximum allowed scattering at minimum :math:`p_{T}` scaled by the forward angle (as
calculated above). Both the check against min :math:`p_{T}` and the check against the
calculated :math:`p_{T}` (discussed further below) take the correlated measurement
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

Following check takes into account estimate particle momentum (smaller scattering
angle is permitted for higher momentum) and pseudorapidity (larger scattering
takes into account amount of the material crosses that takes depends on the angle).

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

The last cut applied in this function is on the transverse impact parameter (or DCA -
distance of closest approach), which is the distance of the perigee of a track from
the interaction region in :math:`mm` of detector radius. It is calculated and cut on
before storing all top SP compatible with both the current middle SP and current
bottom SP. The cut is calculated in the :math:`u/v` plane using the coefficients
:math:`A` and :math:`B`, and the radius of the SP in the middle layer:

$$
d_0 \\leq \\left| \\left( A - B \\cdot r_M \\right) \\cdot r_M \\right|
$$

The bottom SP and middle SP as well as the collection of top SP is passed to
`SeedFilter::filterSeeds_2SpFixed`, whose collected output for the current middle
SP with all compatible bottom SP and top SP is then passed to
`SeedFilter::filterSeeds_1SpFixed`.

SeedFilter::filterSeeds_2SpFixed
--------------------------------

This function assigns a weight (which should correspond to the likelihood that
a seed is good) to all seeds and applies detector specific section of seeds based on weights.
The weight is a “soft cut”, in that it is only
used to discard tracks if many seeds are created for the same middle SP in
`SeedFilter::filterSeeds_1SpFixed`. This process is important to improving computational
performance and the quality of the final track collections by rejecting lower-quality seeds.

The weight can be influenced by:

#. The transverse (:math:`d_{0}`) and longitudinal (:math:`z_{0}`) impact parameters (the higher the distance the smaller the weight)
#. The number of seeds which may belong to the same particle track (:math:`N_{t}`)
#. Optional detector specific cuts.

The transverse impact parameter is multiplied by the configured factor and subtracted from
the weight, as seeds with higher impact parameters are assumed to be less
likely to stem from a particle than another seed using the same middle SP with
smaller impact parameters. The longitudinal impact parameter is subtracted from
the weight if configured.

The number of seeds only differing in top SP which have similar helix radius
and the same sign (i.e. the same charge) is used to increase the weight, as it
means that more than three SPs that may be from the same particle have
been found. The SPs must have a minimum distance in detector radius,
such that SPs from the same layer cannot be counted towards the
increased weight. The number of found compatible seeds is multiplied by a
configured factor and added to the weight.

The optional detector specific cuts can use the weight and the
three SP to apply a hard cut or change the weight of a seed.

The `filterSeeds_2SpFixed` function also includes a fully configurable seed confirmation step that, when enabled
(`seedConfirmation=True`), classifies higher quality seeds as "quality confined" seeds if
they fall within a predefined range of parameters (:math:`d_{0}`, :math:`z_{0}` and :math:`N_{t}`) that also
depends on the region of the detector (i.e., forward or central region). If the seed is not
classified as "quality confined" seed, it will only be accepted if its weight is greater
than a certain threshold and no other high quality seed has been found.

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
The seed is kept only if its weight is greater or equal than the weight of at least one of
its SP components.


Footnotes
---------

.. [#f1] approximately, this is one of the reasons the algorithm breaks down for low energy particles.

