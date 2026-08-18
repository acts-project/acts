@defgroup seeding Seeding
@ingroup pattern_recog
@brief Algorithms and data structures for track seeding

See @ref sec_seeding for the conceptual motivation behind seeding. This page
documents the concrete seeding implementation in ACTS.

A good seeding algorithm has the following properties:

- it finds at least one seed for each particle that should be found.
- it doesn't find many seeds which do NOT correspond to particles.
- it doesn't find many seeds per particle.

## ACTS Implementation

The triplet-based seeding implementation in `Core/include/Acts/Seeding2/` was
written with a focus on parallelism and maintainability, and is as
detector-agnostic as possible, only assuming a (near) homogeneous magnetic
field with particles originating from the central detector region. Cuts are
configurable and can be plugged in as an algorithm called by the seeding. The
seeding works on measurements or "SpacePoints" (SP), which need to provide
@f$(x,y,z)@f$ coordinates with the @f$z@f$ axis along the magnetic field. For
the seeding algorithm to function, the particle's point of origin has to have
a radius smaller than the radius of the detector layer closest to the
interaction region -- in other words, this seeding algorithm is not suitable
for secondary particles originating far from the interaction region.

@anchor fig_seeding_3d

![Sketch of the detector with 3 layers. The interaction region is supposed to be located along the z-axis and has a size significantly smaller than the radius of the innermost detector layer.](seeding/3Dcoordinates.svg) {width=550px}

> [!note]
> The seeding algorithm breaks down for particles whose helix diameter is
> smaller than the detector radius up to which seeds are created. This is due
> to ordering assumptions of SP locations, as well as approximations that
> become inaccurate for lower energy particles.

### SP Grid and Groups Formation

The SPs in each detector layer are projected onto a rectangular grid of
configurable granularity. The search for seeds starts by selecting an SP in
the middle detector layer, then matching SPs are searched for in the inner
and outer layers. Grouping the SPs in this grid limits the search to
neighbouring grid cells, significantly improving algorithm performance (see
@ref fig_triplets_formation "the figure below"). The number of neighboring
bins used in the SP search can be defined separately for the bottom and top
layer SPs in the @f$z@f$ and @f$\phi@f$ directions.

@anchor fig_triplets_formation

![Representation of the search for triplet combinations in the (r, z) plane. The bins used in the search are represented in different colours.](seeding/tripletsFormation.svg) {width=450px}

### The Seed Finder

@ref Acts::TripletSeeder::createSeedsFromGroup (and its multi-group overload
@ref Acts::TripletSeeder::createSeedsFromGroups) orchestrates seed creation
from bottom, middle and top space point ranges constructed from detector
layers of increasing radii. It combines three collaborators:

- @ref Acts::DoubletSeedFinder builds compatible bottom-middle and
  middle-top *doublets* for a given middle SP by applying cuts that can be
  tested with two SPs only (pseudorapidity, origin along the @f$z@f$-axis,
  distance in @f$r@f$ between SPs, compatibility with the interaction
  point).
- @ref Acts::TripletSeedFinder combines compatible bottom and top doublets
  sharing the same middle SP into *triplet* candidates, confronting each
  with the helix hypothesis described below.
- @ref Acts::ITripletSeedFilter (see @ref Acts::BroadTripletSeedFilter for
  the default implementation) ranks and selects the resulting triplet
  candidates, as described in @ref sec_seed_filter "the section below".

In order to perform calculations only once, the circle calculation used for
the helix hypothesis is spread out over the doublet- and triplet-forming
steps rather than repeated for every combination.

@anchor fig_xy_coordinates

![The x-y projection of the detector with the charged particle helical track originating from the centre of the detector. Signals left by the passage of the track through the detector layers are marked with green crosses.](seeding/x-yCoordinates.svg) {width=400px}

From the helix circle (see @ref fig_xy_coordinates "the figure above"),
particle energy and impact parameters can be estimated. To calculate the
helix circle in the @f$x/y@f$ plane, the @f$x,y@f$ coordinates are
transformed into a @f$u/v@f$ plane, to calculate the circle with a linear
equation instead of a quadratic one, for speed. The conformal transformation
is given by:

@f[
u = \frac{x}{x^2+y^2}, \quad \quad v = \frac{y}{x^2+y^2} ,
@f]

where the circle containing the three SPs is transformed into a line with
equation @f$v = Au + B@f$. The angular coefficient @f$A@f$ can be evaluated
by the slope of the linear function between the top and bottom layer SPs,
after transforming their coordinates from @f$x/y@f$ to @f$u/v@f$ using the
previous equations:

@f[
A = \frac{v_t-v_b}{u_t-u_b} ,
@f]

Then, @f$B@f$ can be obtained by inserting @f$A@f$ into the linear equation
for the bottom layer SP:

@f[
v_b = Au_b + B \rightarrow B = v_b - Au_B .
@f]

Inserting the coefficients into the circle equation and assuming the circle
goes through the origin, gives:

@f[
(2R)^2 = \frac{A^2+1}{B^2} .
@f]

This means a cut on the estimate of the minimum helix diameter
`minHelixDiameter2` can be applied without the extra overhead of conversions
or computationally complex calculations. The seed is accepted if

@f[
\frac{A^2+1}{B^2} > (2 R^{min})^2 = \left ( \frac{2 \cdot p_T^{min}}{300 \cdot B_z} \right)^2 \equiv \textnormal{minHelixDiameter2} ,
@f]

where @f$B_z@f$ is the magnetic field.

@anchor fig_rz_coordinates

![The r-z projection of the detector with the same charged particle track. The track is depicted with the same colours as in the previous figure.](seeding/r-zCoordinates.svg) {width=500px}

The track is not an ideal helix: at each detector layer (or any other
material) scattering may occur, making the helix approximate. The algorithm
checks whether the triplet forms a nearly straight line in the @f$r/z@f$
plane (see @ref fig_rz_coordinates "the figure above"), as the particle path
in the @f$r/z@f$ plane is unaffected by the magnetic field. This is split
into two parts. The first test occurs before the calculation of the helix
circle: the deviation from a straight line is compared to the maximum
allowed scattering at minimum @f$p_T@f$, scaled by the forward angle:

@f[
\left (\cot \theta_b - \cot \theta_t \right )^2 < \sigma^2_{p_T^{min}} + \sigma_f^2,
@f]

This check takes into account the squared uncertainty in the difference
between slopes (@f$\sigma_f^2@f$) and the scattering term
`scatteringInRegion2` (@f$\sigma^2_{p_T^{min}}@f$), which is calculated from
`sigmaScattering` (the configurable number of sigmas of scattering angle to
consider) and `maxScatteringAngle2`, evaluated from the Lynch & Dahl
correction (G.R. Lynch and O.I. Dahl, Nucl. Instrum. Methods B58, Phys. Rev.
D 98 6 (1991) 2) of the Highland equation (@cite ParticleDataGroup:2018ovx)
assuming the lowest allowed @f$p_T@f$:

@f[
\Theta_0^{min} = \frac{13.6 \text{MeV}}{p_T^{min}} \cdot \sqrt{\frac{z_q^2 L}{\beta^2 L_0}} \left(1+0.038 \ln \left(\frac{z_q^2 L}{\beta^2 L_0}\right) \right),
@f]

The second part checks against the multiple scattering term
`p2scatterSigma` (@f$\sigma^2_{p_T^{estimated}}@f$), assuming the seed
@f$p_T@f$ estimation instead of the minimum allowed @f$p_T@f$. This term is
inversely proportional to the momentum and accounts for the curvature of the
seed -- a smaller scattering angle is permitted for higher momentum.

@f[
\left (\cot \theta_b - \cot \theta_t \right ) ^2 < \sigma^2_{p_T^{estimated}} + \sigma_f^2,
@f]

The last cut applied in this function is on the transverse impact parameter
(or DCA -- distance of closest approach), the distance of the perigee of a
track from the interaction region in mm of detector radius. It is calculated
and cut on before storing all top SPs compatible with both the current
middle SP and current bottom SP.

@anchor fig_impact_parameter

![Helix representation in the x/y reference frame with the central space-point (SP_m) at the origin.](seeding/impactParameter.svg) {width=400px}

Assuming the middle layer SP is at the origin of the @f$x/y@f$ frame, as in
@ref fig_impact_parameter "the figure above", the distance between the
centre of the helix and the interaction point (IP) is given by

@f[
(x_0 + r_m)^2 + y_0^2 = (R + d_0)^2 \quad  \xrightarrow{R^2 = x_0^2 + y_0^2} \quad \frac{d_0^2}{R^2} + 2 \frac{d_0}{R} = \frac{2 x_0 r_m + r_m^2}{R^2} .
@f]

Considering that @f$d_0 \ll R@f$ (the term proportional to @f$d_0^2@f$ can be
neglected) and using the @f$u/v@f$ line equation calculated previously, the
cut can now be estimated using a linear function in the @f$u/v@f$ plane
instead of a quartic function:

@f[
d_0 \leq \left| \left( A - B \cdot r_M \right) \cdot r_M \right|
@f]

### The Seed Filter {#sec_seed_filter}

After creating the potential seeds, a seed filter procedure compares the
seeds with other SPs compatible with the seed curvature. This process ranks
the potential seeds based on certain quality criteria, and selects the ones
more likely to produce high-quality tracks. The @ref Acts::ITripletSeedFilter
interface (implemented by @ref Acts::BroadTripletSeedFilter for the default
configuration) is divided into two functions, @ref
Acts::ITripletSeedFilter::filterTripletTopCandidates and @ref
Acts::ITripletSeedFilter::filterTripletsMiddleFixed.

The first function compares the middle and bottom layer SPs of the seeds to
other top layer SPs; seeds only differing in the top SP are compatible if
they have a similar helix radius with the same sign (i.e. the same charge).
The SPs must have a minimum distance in detector radius, so that SPs from the
same layer cannot be considered compatible. The second function iterates
over seeds with only a common middle layer SP and selects the higher quality
combinations.

@ref Acts::ITripletSeedFilter::filterTripletTopCandidates assigns a weight
(corresponding to the likelihood that a seed is good) to all seeds and
applies detector-specific selection of seeds based on weight. The weight is
a "soft cut", meaning it is only used to discard tracks if many seeds are
created for the same middle SP. This is important both for computational
performance and for the quality of the final track collections, by
rejecting lower-quality seeds.

The weight is calculated by:

@f[
w = (c_1 \cdot N_{t} - c_2 \cdot d_0 - c_3 |z_0| ) + \textnormal{detector specific cuts}.
@f]

The transverse (@f$d_0@f$) and longitudinal (@f$z_0@f$) impact parameters are
multiplied by a configured factor and subtracted from the weight, since seeds
with higher impact parameters are assumed to be less likely to stem from a
particle than another seed using the same middle SP with smaller impact
parameters. The number of compatible seeds (@f$N_t@f$) increases the weight,
as more measurements lead to higher quality tracks. Finally, the weight can
also be affected by optional detector-specific cuts.

@ref Acts::ITripletSeedFilter::filterTripletTopCandidates also includes a
configurable @ref Acts::SeedConfirmationRangeConfig seed confirmation step
that, when enabled, classifies higher quality seeds as "quality confirmed"
seeds if they fall within a predefined range of parameters (@f$d_0@f$,
@f$z_0@f$ and @f$N_t@f$) that also depends on the region of the detector
(i.e., forward or central). If a seed is not classified as "quality
confirmed", it is only accepted if its weight is greater than a certain
threshold and no other high-quality seed has been found.

The seed confirmation also sets a limit on the number of seeds produced for
each middle SP, retaining only the higher quality seeds. If this limit is
exceeded, the algorithm checks whether there is any low-quality seed in the
seed container of this middle SP that can be removed.

@ref Acts::ITripletSeedFilter::filterTripletsMiddleFixed allows the
detector-specific cuts to filter based on all seeds with a common middle SP,
and limits the number of seeds per middle SP to the configured limit. It
sorts the seeds by weight and, to achieve a well-defined ordering in the rare
case where weights are equal, sorts them by location. The ordering by
location is only done to make sure reimplementations (such as the GPU code)
are comparable and return bitwise exactly the same result.
