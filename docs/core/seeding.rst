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
   :width: 550
   
   Sketch of the detector with 3 layers. The interaction region is supposed to be located along the z axis and have size significantly smaller than the radius of the innermost detector layer.

The SPs in each detector layer are projected on a rectangular grid of configurable
granularity. The search for seed starts from selecting SP in the middle detector 
layer. Then matching SPs are searched in the inner an outer layers. Grouping of 
the SPs in the aforementioned grid allows to limit the search to neighbouring grid
cells thus improving significantly algorithm performance. The number of neighboring 
bins used in the SP search for each :math:`(\phi, z)` bin can be defined
in the :math:`z` direction with the vectors `zBinNeighborsTop` and `zBinNeighborsBottom`
separately for bottom and top SPs, and in the :math:`\phi` direction with `numPhiNeighbors`.

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
two SP only (pseudorapidity, origin along :math:`z`-axis, distance in `r` between SP,
compatibility with interaction point).

If both compatible bottom and top SP have been found, test each bottom SP,
middle SP, top SP triplet combination in a triple nested loop. A major part of
this is the calculation of the helix circle. In order to perform calculations
only once, the circle calculation is spread out over the three loops.

.. figure:: ../figures/seeding/x-yCoordinates.svg
   :name: xy
   :align: center
   :width: 400
   
   The x-y projection of the detector with the charged particle helical track originating from the centre of the detector. Signals left by passage of the track through the detector layers are marked with green crosses.

From the helix circle, particle energy and impact parameters can be estimated.
To calculate the helix circle in the :math:`x/y` plane, the x,y coordinates are
transformed into a :math:`u/v` plane in order to calculate the circle with a linear equation
instead of a quadratic equation for speed. The conformal transformation is given by:
$$
u = \\frac{x}{x^2+y^2}, \\quad \\quad v = \\frac{y}{x^2+y^2} ,
$$
where the circle containing the three SPs are transformed into a line with equation :math:`v = Au + B`.
The angular coefficient :math:`A` can be evaluated by the slope of the linear function between the top and bottom layer SPs, after transforming the coordinates of these SPs from :math:`x/y` to :math:`u/v` using the previous equations:
$$
A = \\frac{v_t-v_b}{u_t-u_b} ,
$$

Then, :math:`B` can be obtained by inserting :math:`A` into the linear equation for the bottom layer SP:
$$
v_b = Au_b + B \\rightarrow B = v_b - Au_B .
$$

Inserting the coefficients in the circle equation and assuming that the circle goes through the origin we obtain:
$$
(2R)^2 = \\frac{A^2+1}{B^2} .
$$

Now we can we can apply a cut on the estimate of the minimum helix diameter (`minHelixDiameter2`) without the extra overhead of conversions or computationally complex calculations. The seed is accepted if
$$
\\frac{A^2+1}{B^2} > (2 R^{min})^2 = \\left ( \\frac{2 \\cdot p_T^{min}}{300 \\cdot B_z} \\right)^2 ,
$$
where :math:`B_z` is the magnetic field.

.. figure:: ../figures/seeding/r-zCoordinates.svg
   :name: rz
   :align: center
   :width: 500
   
   The r-z projection of the detector with the same charged particle track. The track is depicted with the same colours as on previous figure.
       
The next cuts check if the triplet forms a nearly straight line
in the :math:`r/z` plane (see :numref:`rz`) as the particle path in the :math:`r/z` plane is
unaffected by the magnetic field [#f1]_. This is split in two parts; the first test occurs before the calculation of the helix
circle. Therefore, the deviation from a straight line is compared to the
maximum allowed scattering at minimum :math:`p_{T}` scaled by the forward angle:

$$
\\left ( \\frac{1}{\\tan \\theta_b} - \\frac{1}{\\tan \\theta_t} \\right )^2 < \\sigma^2_{p_T^{min}} + \\sigma_f^2,
$$

The scattering term (`scatteringInRegion2` :math:`\equiv \sigma^2_{p_T^{min}}`) is calculated from
`sigmaScattering`, the configurable number of sigmas of scattering angle
to be considered, and `maxScatteringAngle2`, which is evaluated from the
Lynch & Dahl correction of the Highland equation assuming the lowest
allowed :math:`p_{T}`. The parameters of the Highland equation are fully configurable.
The calculation is also spread over the nested loops to avoid redoing calculations.

Following check takes into account estimate particle momentum (smaller scattering
angle is permitted for higher momentum) and pseudorapidity (larger scattering
takes into account amount of the material crosses that takes depends on the angle):

$$
\\left ( \\frac{1}{\\tan \\theta_b} - \\frac{1}{\\tan \\theta_t} \\right ) ^2 < \\sigma^2_{p_T^{estimated}} + \\sigma_f^2,
$$

Both the check against min :math:`p_{T}` and the check against the
calculated :math:`p_{T}` take into account the squared uncertainty in
the difference between slopes (`error2` :math:`\equiv \sigma_f^2`).
By assuming Gaussian error propagation, we can add the two errors
if they are uncorrelated (which is fair for scattering and measurement uncertainties).

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
The weight is a “soft cut”, which means that it is only
used to discard tracks if many seeds are created for the same middle SP in
`SeedFilter::filterSeeds_1SpFixed`. This process is important to improving computational
performance and the quality of the final track collections by rejecting lower-quality seeds.

The weight can be influenced by:

#. The transverse (:math:`d_{0}`) and longitudinal (:math:`z_{0}`) impact parameters (the higher the distance the smaller the weight).
#. The number of seeds which may belong to the same particle track (:math:`N_{t}`).
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

