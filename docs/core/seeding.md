# Track Seeding

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

* It finds at least one seed for each particle that should be found
* It doesn’t find many seeds which do NOT correspond to particles
* It doesn’t find many seeds per particle

The most typical way to create seeds is to combine measurements. In a homogeneous magnetic field, 3 measurements perfectly describe the helical path of a charged particle. One such triplet of measurements would then constitute a seed and define in close bounds where the tracking needs to look for additional measurements to create a track spanning the whole detector. The difficulty is in choosing the correct measurements, as a helix can be fitted through any 3 measurements in a collision event with potentially tens of thousands of measurements. Therefore, a number of constraints or “cuts” are defined to reduce the number of candidates. Cuts may define where particles originate or the range of energy of particles to be found or otherwise restrict the combination of measurements for seed creation.

## Acts Implementation

The seeding implementation in `Core/include/Acts/Seeding/` was written with a focus on parallelism and maintainability and as detector agnostic as possible, only assuming a (near) homogeneous magnetic field with particles originating from the central detector region. Cuts are configurable and can be plugged in as algorithm which is called by the seeding. The seeding works on measurements or “SpacePoints” (SP), which need to provide $(x,y,z)$ coordinates with the $z$ axis being along the magnetic field, and $x$ and $y$. For the seeding algorithm to function the particle point of origin has to have radius smaller compared to the radius of the detector layer closest to the interaction region. In other words, this seeding algorithm is not suitable for secondary particles originating far from interaction region.

:::{figure} ../figures/seeding/3Dcoordinates.svg
:width: 550px
:align: center
Sketch of the detector with 3 layers. The interaction region is supposed to be located along the z axis and have size significantly smaller than the radius of the innermost detector layer.
:::

:::{attention}
Note that the seeding algorithm breaks down for particles with a particle
track whose helix diameter is smaller than the detector radius until which
seeds are to be created. This is due to ordering assumptions of SP
locations as well as due to approximations which become inaccurate for
lower energy particles.
:::

### SP Grid and Groups Formation

The SPs in each detector layer are projected on a rectangular grid of configurable
granularity. The search for seed starts from selecting SP in the middle detector 
layer. Then matching SPs are searched in the inner an outer layers. Grouping of 
the SPs in the aforementioned grid allows to limit the search to neighbouring grid
cells thus improving significantly algorithm performance (see {numref}`tripletsFormation`). The number of neighboring 
bins used in the SP search can be defined separately for bottom and top layer SPs
in the $z$ and $\phi$ directions.

(tripletsFormation)=
:::{figure} ../figures/seeding/tripletsFormation.svg
:width: 450px
:align: center
Representation of the search for triplet combinations in the $(r, z)$ plane. The bins used in the search are represented in different colours.
:::

### The Seed Finder

The The {class}`SeedFilter::createSeedsForGroup` function receives three iterators 
over SPs constructed from detector layers of increasing radii. The seedfinder will 
then attempt to create seeds, with each seed containing exactly one SP returned by 
each of the three iterators. It starts by iterating over SPs in the middle layer
(2nd iterator), and within this loop separately iterates once over bottom SP 
and once over top SP. Within each of the nested loops, SP pairs are tested for
compatibility by applying a set of configurable cuts that can be tested with
two SP only (pseudorapidity, origin along $z$-axis, distance in $r$ between SP,
compatibility with interaction point).

For all pairs passing the selection the triplets of bottom - middle - top SPs are formed.
Each triplet is then confronted with the helix hypothesis. In order to perform calculations
only once, the circle calculation is spread out over the three loops.

(x-yCoordinates)=
:::{figure} ../figures/seeding/x-yCoordinates.svg
:width: 400px
:align: center
The x-y projection of the detector with the charged particle helical track originating from the centre of the detector. Signals left by passage of the track through the detector layers are marked with green crosses.
:::

From the helix circle (see {numref}`x-yCoordinates`), particle energy and impact parameters can be estimated.
To calculate the helix circle in the $x/y$ plane, the $x,y$ coordinates are
transformed into a $u/v$ plane in order to calculate the circle with a linear equation
instead of a quadratic equation for speed. The conformal transformation is given by:

\begin{equation*}
u = \frac{x}{x^2+y^2}, \quad \quad v = \frac{y}{x^2+y^2} ,
\end{equation*}

where the circle containing the three SPs are transformed into a line with equation $v = Au + B$.
The angular coefficient $A$ can be evaluated by the slope of the linear function between the top and bottom layer SPs, after transforming the coordinates of these SPs from $x/y$ to $u/v$ using the previous equations:

\begin{equation*}
A = \frac{v_t-v_b}{u_t-u_b} ,
\end{equation*}

Then, $B$ can be obtained by inserting $A$ into the linear equation for the bottom layer SP:
\begin{equation*}
v_b = Au_b + B \rightarrow B = v_b - Au_B .
\end{equation*}

Inserting the coefficients in the circle equation and assuming that the circle goes through the origin we obtain:

\begin{equation*}
(2R)^2 = \frac{A^2+1}{B^2} .
\end{equation*}

Now we can we can apply a cut on the estimate of the minimum helix diameter `minHelixDiameter2` without the extra overhead of conversions or computationally complex calculations. The seed is accepted if

\begin{equation*}
\frac{A^2+1}{B^2} > (2 R^{min})^2 = \left ( \frac{2 \cdot p_T^{min}}{300 \cdot B_z} \right)^2 \equiv \textnormal{minHelixDiameter2} ,
\end{equation*}

where $B_z$ is the magnetic field.

(r-zCoordinates)=
:::{figure} ../figures/seeding/r-zCoordinates.svg
:width: 500px
:align: center
The r-z projection of the detector with the same charged particle track. The track is depicted with the same colours as on previous figure.
:::
       
The track is not an ideal helix. At each detector layer (or any other material)
scattering may occur making the helix approximate. 
The algorithm will check if the triplet forms a nearly straight line
in the $r/z$ plane (see {numref}`r-zCoordinates`) as the particle path in the $r/z$ plane is
unaffected by the magnetic field. This is split in two parts; the first test occurs before the calculation of the helix
circle. Therefore, the deviation from a straight line is compared to the
maximum allowed scattering at minimum $p_T$ scaled by the forward angle:

\begin{equation*}
\left (\cot \theta_b - \cot \theta_t \right )^2 < \sigma^2_{p_T^{min}} + \sigma_f^2,
\end{equation*}

This check takes into account the squared uncertainty in the difference between slopes ($\sigma_f^2$) and the
scattering term `scatteringInRegion2`  ($\sigma^2_{p_T^{min}}$), which is calculated from
`sigmaScattering`, the configurable number of sigmas of scattering angle
to be considered, and `maxScatteringAngle2`, which is evaluated from the
Lynch & Dahl correction[^Tanabashi:2018]$^,$[^Lynch:1991] of the Highland equation assuming the lowest
allowed $p_T$:

\begin{equation*}
\Theta_0^{min} = \frac{13.6 \text{MeV}}{p_T^{min}} \cdot \sqrt{\frac{z_q^2 L}{\beta^2 L_0}} \left(1+0.038 \ln \left(\frac{z_q^2 L}{\beta^2 L_0}\right) \right),
\end{equation*}

The second part check against the multiple scattering term `p2scatterSigma` ($\sigma^2_{p_T^{estimated}}$) assuming the seed $p_T$ estimation, instead of the minimum allowed $p_T$.
This term is inversely proportinal to the momentum; smaller scattering angle is permitted for higher momentum, while larger scattering
takes into account amount of the curvature of the seed:

\begin{equation*}
\left (\cot \theta_b - \cot \theta_t \right ) ^2 < \sigma^2_{p_T^{estimated}} + \sigma_f^2,
\end{equation*}

The last cut applied in this function is on the transverse impact parameter (or DCA -
distance of closest approach), which is the distance of the perigee of a track from
the interaction region in $mm$ of detector radius. It is calculated and cut on
before storing all top SP compatible with both the current middle SP and current
bottom SP. 

(impactParameter)=
:::{figure} ../figures/seeding/impactParameter.svg
:width: 400px
:align: center
oi
:::

Assuming the middle layer SP in the origin of the $x/y$ frame, as in {numref}`impactParameter`. 
The distance between the center of the helix and interaction point (IP) is given by
\begin{equation*}
(x_0 + r_m)^2 + y_0^2 = (R + d_0)^2 \quad  \xrightarrow{R^2 = x_0^2 + y_0^2} \quad \frac{d_0^2}{R^2} + 2 \frac{d_0}{R} = \frac{2 x_0 r_m + r_m^2}{R^2} .
\end{equation*}

Considering that $d_0 << R$ (we can neglect the term proportional to $d_0^2$) and using the $u/v$ line equation calculated previously, 
the cut can now be estimated using a linear function in the $u/v$ plane instead of a quatic function:

\begin{equation*}
d_0 \leq \left| \left( A - B \cdot r_M \right) \cdot r_M \right|
\end{equation*}

### The Seed Filter

The bottom SP and middle SP as well as the collection of top SP is passed to
{class}`SeedFilter::filterSeeds_2SpFixed`, whose collected output for the current middle
SP with all compatible bottom SP and top SP is then passed to
{class}`SeedFilter::filterSeeds_1SpFixed`.

#### SeedFilter::filterSeeds_2SpFixed

This function assigns a weight (which should correspond to the likelihood that
a seed is good) to all seeds and applies detector specific section of seeds based on weights.
The weight is a “soft cut”, which means that it is only
used to discard tracks if many seeds are created for the same middle SP in
{class}`SeedFilter::filterSeeds_1SpFixed`. This process is important to improving computational
performance and the quality of the final track collections by rejecting lower-quality seeds.

The weight can be influenced by:

* The transverse ($d_0$) and longitudinal ($z_0$) impact parameters (the higher the distance the smaller the weight).
* The number of seeds which may belong to the same particle track ($N_t$).
* Optional detector specific cuts.

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

The {class}`SeedFilter::filterSeeds_2SpFixed` function also includes a fully configurable seed confirmation step that, when enabled,
classifies higher quality seeds as "quality confined" seeds if
they fall within a predefined range of parameters ($d_0$, $z_0$ and $N_t$) that also
depends on the region of the detector (i.e., forward or central region). If the seed is not
classified as "quality confined" seed, it will only be accepted if its weight is greater
than a certain threshold and no other high quality seed has been found.

The seed confirmation also sets a limit on the number of seeds produced for each middle SP,
which retains only the higher quality seeds. If this limit is exceeded, the algorithm
checks if there is any low-quality seed in the seed container of this middle SP that can be removed.

#### SeedFilter::filterSeeds_1SpFixed

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


[^Tanabashi:2018]: M. Tanabashi et al. (Particle Data Group), Passage of Particles Through Matter, Phys. Rev. D 98 0300001 (2018) 2.

[^Lynch:1991]: G.R. Lynch and O.I Dahl, Nucl. Instrum. Methods B58, Phys. Rev. D 98 6 (1991) 2.
