@page tracking Tracking in a nutshell

Track reconstruction is the process to recover the properties of a charged
particle from a set of measurements caused by the interaction with some form of
sensitive detector. The goal is to find which measurements are likely to have
been caused by which particle, to group them accordingly, and to estimate the
associated trajectory. Such charged particle trajectories form the basic input
to the majority of higher-level reconstruction procedures in many cases.

![Illustration of a track reconstruction chain starting from space points to fully formed tracks.](tracking/tracking.svg) {width=60%}

This section provides a high-level view of a track reconstruction chain, and is
largely based on @cite gessinger_befurt_2021_m23d2-xsq75. It gives an overview of the basic
building blocks conceptually, and also connects these building blocks with the
concrete implementations in the core ACTS library, where available.

# Charged particle detection {#charged-particle-detection}

The first step in the chain to reconstruct charged particles is their
detection using sensitive elements. Charged particle detection can be
achieved in a variety of ways with very different technologies. It can be
achieved by measuring the interaction of charged particles with matter. When
this occurs, the interacting particle typically ionizes the surrounding
material. Particle detectors make use of this fact by converting the
resulting charge into a measurable signal in various ways.


@anchor segmentation

![Illustration of a one-dimensional (a) and a two-dimensional segmentation (b) of a silicon sensor.](tracking/segmentation.svg) {width=70%}

A very common electronic detection approach is the use of semiconducting
particle detectors, often made of silicon. When a charged particle traverses
such a sensor, it ionizes the material in the depletion zone, caused by the
interface of two different semiconducting materials. The result are pairs of
opposite charges. These charge pairs are separated by an electric field and
drift toward the electrodes. At this point, an electric signal is created
which can be amplified and read out. By means of segmentation, the measured
signal can be associated with a location on the sensor. Silicon sensors are
usually segmented in one dimension (*strips*) or in two dimensions
(*pixels*) (see Figure @ref segmentation).

# Track parametrization {#track-parametrization}

To express the properties of a particle's trajectory, a choice of parameters
has to be made. The parameters need to express all the relevant quantities of
interest. In the presence of a magnetic field, which affects charged
trajectories, the global position and momentum, as well as the charge are
needed to fully specify the particle properties. In addition, a time parameter
can be included. Apart from the global reference frame, track quantities often
need to be represented with respect to a surface. This can be achieved with a
parametrization like

@f[
  \vec x = \left(l_0, l_1, \phi, \theta, q/p, t\right)^T
@f]

although other parameter conventions exist as well.
Figure @ref parameters illustrates this choice of parameters. @f$l_0@f$, @f$l_1@f$ are
the local coordinates of the corresponding surface, @f$\phi \in [-\pi,\pi)@f$ and
@f$\theta \in [0,\pi]@f$ are the angles in the transverse and longitudinal
direction of the global frame, expressed with respect to the current location
along the trajectory, as indicated in Figure @ref parameters (b). @f$\theta@f$ is
the polar angle relative to the positive @f$z@f$-axis, and @f$\phi@f$ is the azimuth
angle in the transverse plane. Finally, @f$q/p@f$ combines the charge of the
particle with the inverse momentum. In Figure @ref parameters (a), the global
momentum vector @f$\vec p@f$ is shown, which can be recovered from the parameters
@f$\vec x@f$ using @f$\phi@f$, @f$\theta@f$ and @f$q/p@f$.

@anchor parameters

![Illustration of the parametrization of a particle track with respect to a two-dimensional surface. (a) shows the local position, global momentum and their corresponding uncertainties. (b) displays the angles $\\phi$ and $\\theta$ in the transverse and longitudinal planes.](tracking/parameters.svg) {width=70%}

@anchor perigee

![Illustration of the perigee parametrization which uses the point of closest approach relative to a reference point. The impact parameter $d_0$, the position $l$ and the momentum vector $\\vec p$ are shown.](tracking/perigee.svg) {width=40%}

Aside from the nominal quantities captured in @f$\vec x@f$, the related
uncertainties and correlations need to be taken into account as well. They
can be expressed as a @f$6\times 6@f$ covariance matrix like

@f[
  C =
  \begin{bmatrix}
   \sigma^2(l_0)& \text{cov}(l_0,l_1) & \text{cov}(l_0, \phi) & \text{cov}(l_0, \theta) & \text{cov}(l_0, q/p) & \text{cov}(l_0, t) \\
   . & \sigma^2(l_1) & \text{cov}(l_1, \phi) & \text{cov}(l_1, \theta) & \text{cov}(l_1, q/p) & \text{cov}(l_1, t) \\
   . & . &  \sigma^2(\phi) & \text{cov}(\phi,\theta) & \text{cov}(\phi, q/p) & \text{cov}(\phi, t) \\
   . & . & . & \sigma^2(\theta) & \text{cov}(\theta, q/p) & \text{cov}(\theta, t) \\
   . & . & . & . & \sigma^2(q/p) & \text{cov}(q/p, t) \\
   . & . & . & . & . & \sigma^2(t)
  \end{bmatrix}
@f]

Here, @f$\text{cov}(X,Y)@f$ is the covariance of variables @f$X@f$ and @f$Y@f$, while
@f$\sigma^2(X)@f$ are the regular variances. As the covariance matrix @f$C@f$ is
symmetric, only the upper right half is shown in the matrix above. The
uncertainties associated with the local position, as well as the momentum
direction are indicated in Figure @ref parameters (a) as an ellipse and a cone
around the momentum vector @f$\vec p@f$, respectively.

# Particle propagation {#particle-propagation}

> [!tip]
> A dedicated description of the ACTS implementation of particle propagation
> can be found @ref propagation "here".

A central part of track reconstruction is the ability to calculate the
trajectory of a charged particle, given its properties at a given point. This
process, called *particle propagation* or *extrapolation*, is used to predict a
particle's properties after it has travelled a certain distance. In many cases,
the projected intersection with various types of surfaces is desired. The
trajectory of a charged particle is governed by the @ref magnetic_field
"magnetic field" through which it travels, as well as any material effects (see
@ref material). In case of a homogeneous magnetic field, and in the absence of
material interaction, the particle follows a helical trajectory. Such a helix
can be calculated purely analytically, although intersections require numerical
methods nevertheless.

Often, magnetic fields are not homogeneous, however. In the presence of such
changing fields, the corresponding differential equations of motions need to be
solved using numerical integration techniques.

## Numerical integration {#numerical-integration}

In ACTS, numerical integration is done using the *Runge-Kutta-Nyström* (RKN) method.
Commonly used in the variant at fourth order, the RKN method describes how to calculate a
solution to an initial value problem that can be formulated generically like


@f[
\frac{dy}{dt} = f(t,y), \qquad y(t_0) = y_0,
@f]


where @f$y_0@f$ refers to the initial value of @f$y@f$ at @f$t_0@f$, and
@f$f(t,y)@f$ is the functional form describing the dynamics. The method then
successively approximates the analytical solution @f$y(t)@f$ in a stepwise
fashion. At each step @f$(t_n, y_n)@f$, the goal is effectively to approximate
the next step @f$y(t_{n+1})@f$. Using a step size @f$h@f$, the algorithm evaluates
the function @f$f@f$ at four points @f$k_{1-4}@f$:


@f[
\begin{aligned}
    k_1 &= f(t_n, y_n) \\
    k_2 &= f\left( t_n + \frac h 2, y_n + h \frac{k_1} 2 \right) \\
    k_3 &= f\left( t_n + \frac h 2, y_n + h \frac{k_2} 2 \right)\\
    k_4 &= f\left( t_n + h, y_n + hk_3 \right).
\end{aligned}
@f]

@anchor rk

![Illustration of the RKN method approximating a first order differential equation. Shown is the calculation of an estimate $y_{n+1}$ at $t_{n+1} = t_n + h$, based on the current step $(t_n,y_n)$. Shown are the four distinct points at which function $y(t)$ is evaluated, and which are blended to form the estimate.](tracking/rk.svg) {width=40%}


Looking at Figure @ref rk, the meaning of these four points in relation
to the step size @f$h@f$ can be understood. @f$k_1@f$ is the derivative at the
current location, @f$k_{2,3}@f$ use @f$k_1@f$ and @f$k_2@f$ respectively to calculate two
envelope derivatives at @f$h/2@f$ and @f$k_4@f$ uses @f$k_3@f$ to make an estimate of the
derivative at @f$h@f$. Combining @f$k_{1-4}@f$, @f$(t_{n+1},y_{n+1})@f$ can be calculated
as the approximation of @f$y(t_{n+1})@f$ like


@f[
\begin{aligned}
    y_{n+1} &= y_n + \frac 1 6 h ( k_1 + 2 k_2 + 2 k_3 + k_4)\\\\
    t_{n+1} &= t_n + h
\end{aligned}
@f]


by effectively averaging the four derivatives. It is apparent that
the step size crucially influences the accuracy of the approximation. A large
step size weakens the approximation, especially if the magnetic field changes
strongly. On the other hand, a too small step size will negatively affect the
execution time of the algorithm.

The Runge-Kutta-Nyström method from above can be adapted to handle second order
differential equations, as is needed for the equations of motion in question,


@f[
\frac{d^2 \vec r}{ds^2} = \frac q p \left( \frac{d\vec r}{ds} \times \vec B (\vec r) \right) = f(s, \vec r, \vec T), \qquad \vec T \equiv \frac{d \vec r}{ds},
@f]



with the global position @f$\vec r@f$, the path element @f$s@f$, the
normalized tangent vector @f$\vec T@f$ and the magnetic field @f$\vec B(\vec r)@f$ at
the global position. A slight modification of @f$k_{1-4}@f$ is also required,
incorporating the first derivative of @f$f(s, \vec r, \vec r')@f$, finally
leading to


@f[
\begin{aligned}
  \vec T_{n+1} &= \vec T_n + \frac h 6 (k_1 + 2k_2 + 2k_3 + k_4) \\\\
  \vec r_{n+1} &= \vec r_n + h \vec T_n + \frac{h^2}{6} (k_1 + k_2 + k_3).
\end{aligned}
@f]


A strategy exists to dynamically adapt the step size according to the magnetic
field strength, with the definition of a target accuracy that the algorithm
tries to achieve. Here, the step size @f$h@f$ will successively be decreased and
the approximation recalculated until the accuracy goal is achieved. Even with
these additional calculations, the approach is still preferable over a
consistently low step size.

## Covariance transport {#covariance-transport}

Aside from the prediction of the track parameters at a given path length, a key
ingredient to many dependent applications are the uncertainties in the form of
the associated covariance matrix @f$C@f$. Conversions between covariance matrices
@f$C^i\to C^f@f$ can generally be achieved like


@f[
C^f = J \cdot C^i \cdot J^T,
@f]


using the Jacobian matrix



@f[
J = \begin{bmatrix}
    \frac{\partial l_0^f}{\partial l_0^i} & \cdots  & \frac{\partial l_0^f}{\partial (t_{0})^i} \\
    \vdots & \ddots & \vdots \\
    \frac{\partial (t_{0})^f}{\partial l_0^i} & \cdots  & \frac{\partial (t_0)^f}{\partial (t_0)^i}
  \end{bmatrix},
@f]



between initial and final parameters @f$\vec x^i@f$ and @f$\vec x^f@f$. The
task therefore becomes calculating the necessary Jacobians to achieve correct
transformation.

One part is the transformation between different coordinate systems, but at the
same location along the trajectory. For this purpose, generic Jacobians can be
calculated between each coordinate system type, and a common coordinate system.
The common coordinate system used for this purpose is the curvilinear frame,
which consists of the global direction angles, and a plane surface located at
the track location, with the normal of the plane surface aligned with the track
momentum. By using Jacobians to the curvilinear frame and the corresponding
inverse matrices, conversions between any two coordinate systems can be
performed.

The second part is the calculation of the evolution of the covariance matrix
during the propagation between surfaces. To this end, a semi-analytical
method which calculates the effective derivatives between two consecutive
RKN steps can be used. By accumulating the Jacobian
matrices calculated for each step, the effective Jacobian between the
starting point and the destination can be obtained.

# Material effects {#material-eff}

Charged particles interact with matter as they pass through it. Since
particle detectors inevitably consist of some form or material, this effect
cannot be completely avoided. By building tracking detectors as light as
possible, and arranging passive components, such as services and support
structures carefully, the material a particle encounters before being
measured can be reduced. Charged particles traversing any form of matter
undergo elastic and inelastic interactions with the atomic structure of the
material, depending on the particle properties.

@anchor multiple_scattering

![Illustration of the effect of multiple scattering on the trajectory of a charged particle passing through a block of material. Entering from the left, it undergoes a series of scattering events, deflecting the trajectory statistically, before exiting on the right.](tracking/multiple_scattering.svg) {width=200px}

In elastic interactions, the particle does not lose a significant amount of
energy, while its trajectory is affected. Figure @ref multiple_scattering shows a
sketch of the way multiple Coulomb scattering affects the direction of a
particle trajectory. In addition, a shift in the transverse plane relative to
the incident direction can occur. As the scattering events occur in
statistically independent directions, the means of both the deflection and
offset tends toward zero as the number of scatters increases. Therefore, in the
numerical particle propagation, this can be accounted for by simply increasing
the uncertainties associated with the direction, depending on the amount of
material encountered.

On the other hand, there are interactions during which the particle loses some
of its energy. Relevant processes here are ionization, as well as
bremsstrahlung for light particles like electrons. For hadronic particles,
hadronic interactions with the nuclei of surrounding material is another
process of interest. In such hadronic interactions, the incoming particle often
disintegrates, and does not propagate further. Since the size of ionization
losses only fluctuates to a small degree for thin layers of material, they can
usually be accounted for by reducing the trajectory energy correspondingly. For
bremsstrahlung, where fluctuations are much larger and the effect cannot be
modelled adequately with a Gaussian distribution, dedicated techniques are
needed (see @ref Acts::GaussianSumFitter).

Two main approaches are implemented in ACTS. The first approximates the
material interaction by using a description that averages the real material
onto thin surfaces across the detector (more on this in
@ref geometry-and-material-modelling). When the propagation encounters such a
surface, it retrieves the material properties, and executes parametrized
modifications to the particle properties and uncertainties. In the second
approach, material effects are continuously incorporated during propagation,
rather than at discrete locations. The latter approach is especially suited for
propagation through volumes of dense material, where the discretization of the
material distribution will not work as well.

# Geometry and material modelling {#geometry-and-material-modelling}

> [!tip]
> A dedicated description of the ACTS implementation of the tracking geometry
> model can be found @ref geometry "here".

A detailed model of the geometry of an experiment is required for tracking. In
many cases, external information is needed to associate a sensitive element
with a position and rotation in the laboratory frame. In case of silicon
sensors, the intrinsic information captured by the sensor is restricted to the
measurement plane. Using a transformation matrix, this local measurement can be
turned into a global one.

Full simulation using tools like Geant4 are frequently used in HEP.
It includes its own geometry
description framework. For the precise simulation of particle interactions
with the detector, this geometry modelling
is highly detailed. Even very small details of the
physical hardware can be crucial, and are often included in the geometry
description. An example for this are readout chips on silicon sensors, or cooling
elements. Figure @ref geometry_detail (a) shows a sketch of such a detailed
geometry description. Shown as an example is a *layer* of silicon
sensors in a barrel configuration. The green rectangles represent the actual
sensitive surfaces, while other elements include cooling, readout and other components.

@anchor geometry_detail

![Sketch of the way a fully detailed simulation geometry (a) models passive elements, in addition to the sensitive elements shown in green. (b) shows a simplified version, where all non-sensitive elements are approximated.](tracking/geometry_detail.svg) {width=90%}

In the majority of cases in track reconstruction, this detailed
geometry is unnecessary. During track reconstruction, the aforementioned
associated information needs to be accessible for measurements, so all
sensitive elements need to be included in some form. Passive elements, on the
other hand, are only required to factor in material interaction effects (see
@ref particle-propagation). Moreover, the fully detailed geometry comes
at the disadvantage of introducing significant overhead during navigation. In
this process, an algorithm attempts to figure out which elements the particle
propagation needs to target, as the trajectory is likely to intersect them.
With a geometry description this precise, the navigation process becomes a
significant performance bottleneck.

@anchor layer_barrel

![Sketch of the way sensitive elements are grouped into layers. Shown is an $xy$-view of a number of sensors, arranged as in e.g. the ATLAS silicon detector barrels. The grouping is based on their mounting radius. The layers are indicated in different colors.](tracking/layer_barrel.svg) {width=40%}

As a compromise between modelling accuracy and performance, ACTS
uses a simplified geometry model. It
focusses on the sensitive elements, which are strictly needed, while passive
elements are discarded from the explicit description and approximated.
Figure @ref geometry_detail (b) shows such a simplified geometry. Here, the
sensitive elements are still shown in green, and other elements are greyed
out, indicating that they are discarded. The sensitive elements are then
grouped into layers, as sketched in Figure @ref layer_barrel. How exactly the
grouping occurs depends on the concrete experiment geometry. In some cases, the layers have the shape
of cylinder surfaces with increasing radii. This example is shown in the
figure in the transverse plane at radii @f$r_{1,2,3}@f$. In the endcaps, where
modules are often arranged on disks, these are used as the layer shape. An
illustration of endcap disk layers can be found in Figure @ref layer_ec, where
six disks are located at six distinct positions in @f$\pm z_{1,2,3}@f$, and
shown in different colors.

@anchor layer_ec

![Sketch of the way sensitive elements are grouped into layers. Shown is a view of a number of sensors, arranged as in e.g. the ATLAS silicon detector endcaps. They are grouped into disks based on their mounting position in $z$. The layers are indicated in different colors.](tracking/layer_ec.svg) {width=80%}

During particle propagation, the navigation makes use of this layer
system. Each layer contains a binned structure, which maps a bin to a set
of sensitive surfaces that overlap with the bin area. This is illustrated in
Figure @ref geo_binning, where the left picture shows the sensitive surface
structure of an exemplary endcap disk. The picture on the right overlays the
binning structure that can be used to enable fast retrieval of compatible
sensitive surfaces. By performing a simple bin lookup, the navigation can
ascertain which sensors it needs to attempt propagation to.


@anchor geo_binning

![Illustration of the binning structure that is used to subdivide layer surfaces. (a) shows two sensor rings of different radii grouped into one disk layer. (b) overlays the binning structure that the navigation queries for compatible surfaces.](tracking/surface_array.svg) {width=80%}

Furthermore, layers are grouped into volumes. Each volume loosely corresponds
to a region of the detector.
Volumes are set up such that their boundary surfaces always touch another
volume. An exception to this is the outermost volume. Each volume's boundary
surfaces store which volume is located on their other side, essentially
forming portals between the volumes. This glueing enables the geometry
navigation between volumes. When the propagation has finished processing a
set of layers, it attempts to target the boundary surfaces. Once a boundary
surface is reached, the active volume is switched, and the next set of layers
is processed.

Care has to be taken to correctly model the passive material, that is
initially discarded with non-sensitive elements. For the material effects to
be correctly taken into account during particle propagation, the material is
projected onto dedicated material surfaces. These material surfaces are
spread across the detector geometry. Each layer is created with two
*approach surfaces* on either side. Their distance can be interpreted as
the thickness of the layer in question. Examples of these approach surfaces
can be found in Figure @ref geometry_detail, at the inner and outer radius.
Approach surfaces, and the boundary surfaces between volumes mentioned before,
are candidates to receive a projection of the surrounding material.
Additional artificial material layers can also be inserted to receive
projected material.

The projection procedure (see @ref material and @ref material_mapping) works
by extrapolating test particles using the fully detailed simulation geometry.
During the extrapolation, the material properties of the geometry are sampled
in small intervals. Subsequently, the same test particle is extrapolated
through the tracking geometry. All material samples are then assigned and
projected onto the closest material surface. Finally, the projection is
averaged. The exact number and placement of the material surfaces has to be
optimized to yield a sufficiently accurate representation of the inactive
material in the detector.

The numerical integration uses these projected material surfaces. Whenever
such a surface is encountered in the propagation, the material properties are
retrieved, and the corresponding modifications to the trajectory are
executed. In case material is supposed to be integrated in a continuous way
(as mentioned in @ref particle-propagation), volumes can also store an
effective volumetric material composition, which is queried by the numerical
integration when needed. As the actual physical location of the detection
hardware can vary over time, possible misalignment of the sensors needs to be
handled correctly.

# Clustering {#clustering}

> [!tip]
> See @ref clustering for information of the implementation of clustering
> in the core library.

The actual track reconstruction procedure itself starts with the conversion of
raw inputs that have been read out from the detector. In case of silicon
detectors, the readout can either be performed in a binary way, only recording
which segments fired, or the amount of charges measured in the segment can be
recorded, e.g. via *time-over-threshold* readout. In all cases, the readout is
attached to an identifier uniquely locating the segment on the corresponding
sensor.

As a next step, these raw readouts need to be *clustered*, in order to
extract an estimate of where particles intersect with the sensor. The general
strategy of clustering algorithms follows the Connected Component Analysis (CCA)
approach, where subsets of segments are successively grouped into clusters.
In case of the Pixel detector, this clustering occurs in two dimensions,
corresponding to the segmentation of its sensors. Here, the CCA can
either consider all eight surrounding pixels as neighboring a central one, or
only consider the four non-diagonal ones, as shown in
Figure @ref clustering_cca. The figure only shows the simplest possible
cluster starting from the central pixel. In reality, the CCA will iteratively
continue from the pixels on the cluster edges.

@anchor clustering_cca

![Illustration of both eight and four cell connectivity.](tracking/cca.svg) {width=60%}

Subsequently, the effective cluster position needs to be estimated. Multiple
factors play a role here. First of all, the average position of the cluster
can be calculated either using only the geometry position of the segments,


@f[
\vec r = \frac{1}{N} \sum_{i=1}^N \vec l_i,
@f]


or be weighted by the charge collected in each segment:


@f[
\vec r = \frac{1}{\sum_{i=1}^N q_i} \sum_{i=1}^N q_i \vec l_i.
@f]


Here, @f$\vec l_i@f$ is the local position of the @f$i@f$-th segment while
@f$q_i@f$ is its charge.

An illustration of the clusterization can be found in Figure @ref clustering_image,
where a pixel sensor is shown to be intersected by a charged particle,
entering on the lower left and exiting on the top right. Three cells shown
with a red frame receive energy from the particle, but the amount is under
the readout threshold. Four other cells receive energy above the threshold
and are read out. The clustering will then group these four cells into a
cluster, and subsequently estimate the cluster position based on the energy
deposited in the cells. In case no charge information is not available
for a given detector, the calculation is purely geometric.


@anchor clustering_image

![Illustration of the clustering of multiple pixels into a cluster, in a three-dimensional view on the left and a projection onto the $xy$-plane on the right. A particle enters the center in the lower left, crosses several segments before exiting the sensor on the top right. The cell colors indicate how far along the trajectory they are encountered.](tracking/clustering.svg) {width=50%}

Another factor that needs to be accounted for is the drift direction of the
created charges. In addition to the collection field of the sensor itself,
the surrounding magnetic field modifies the drift direction by the
*Lorentz-angle* @f$\theta_\text{L}@f$. Depending on the field strength, this
additional angle can cause segments to be activated that would otherwise not
be geometrically within reach of the charges. Other effects, such as the fact
that the modules are not perfectly flat, as the geometry description assumes,
or cross-talk between readout channels, also play a role at this stage.

In the presence of high event activity, particle intersections on single
sensors can be close enough to one another to result in clusters that are not
clearly separated from each other. This circumstance can be somewhat
mitigated by allowing tracks to share clusters with other particles, which
comes at the price of allowing duplicated tracks to some extent.
Additionally, merged clusters typically feature worse position resolution,
which manifests itself since it negatively affects the final fit of the
track.

# Space point formation {#space-point-formation}

The basic input to most forms of pattern recognition algorithms for tracking
are space points, which need to be assembled from the raw measurements. To this
end, the raw measurements are combined with information provided by the
geometry description, such as the location and rotation of the sensors. In this
way, the locations, which are restricted to be local to the sensor surfaces
intrinsically, can be converted into three dimensional points in space.  See
@ref sp_formation for a description of the implementation of
space point formation in the core library.

The @ref fig_sensor "figure" below shows an illustration of the information that is consumed for
a pixel measurement. Shown are three clusters on a sensor, which are caused by
three tracks intersecting it. The corresponding cluster positions are indicated
as well, and can be converted to global positions using the inverse of the
global-to-local transformation matrix, that is provided by the geometry
description.

@anchor fig_sensor

![Illustration of a pixel sensor and its local coordinate system in relation to the global laboratory frame. A transformation allows conversion between the two systems. Shown are three tracks intersecting the sensor, alongside clusters that they produce.](tracking/sp_l2g.svg) {width=50%}

In strip detectors, on the other hand, only a single
dimension is segmented, and an individual measurement is therefore only
constrained in one direction on the surface. However, usually the
strip modules are mounted in pairs, with a stereo angle rotation
between the pairs. To form global space points, measurements from both
sensors of a pair need to be combined.
Due to the stereo angle, a two dimensional
location on the orthogonal projection plane relative to the two parallel
pairs can be found. Using the global transformations of the pair, the
combined measurement location can be converted to global coordinates. If
multiple measurements are located on a stereo pair of strip sensors, there
exists an ambiguity on how to combine strips to form space points, which has to be resolved.


# Seeding {#sec_seeding}

The next step after space point formation is pattern recognition, which be
implemented in various ways.  Global methods exist which attempt to cluster
space points, such as conformal mapping. In this approach, the space points are
transformed into a feature parameter space that reveals patterns for hits
belonging to the same track. In the specific example of a Hough transform, a
parameter space @f$\left(\phi, q/p_\mathrm{T}\right)@f$ is used. As a result, each
space point is effectively transformed into a line, as a series of combinations
of these parameters would lead to the same space point. The lines from a set of
space points of a single track will intersect in one common area. Such an
intersection can be used to identify which space points originate from the same
track. However, this task grows in complexity as detector activity increases
and is susceptible to material effects. See @ref seeding for a description
of the seeding implementation in the core library.

Another group of approaches is the one of seeding and track following. These
algorithms differ from the global ones in that they evaluate individual
combinations of space points, and successively explore the events. One
algorithm from this group is the cellular automaton that iteratively forms
chains of space points going from one layer to the next.

The main approach in ACTS is an algorithm that operates on coarse
subdivisions of the detector is used. This seeding algorithm attempts to find
triplets of space points from increasing radii which are likely to belong to
the same track. It achieves this by iterating the combinatorial triplets and
successively filtering them. Filtering is performed based on the momentum and
impact parameters, which the algorithm attempts to estimate for each triplet.

Under the assumption of a homogeneous magnetic field along the @f$z@f$-axis,
charged particles should follow helical trajectories. In the transverse plane,
the motion is circular, while it is a straight line in the @f$rz@f$-plane.  The
transverse impact parameter and momentum can be estimated from the radius of
the circle in the transverse plane like


@f[
d_0 = \sqrt{c_x^2 + c_y^2} - \rho,
@f]


with the circle center @f$(c_x, c_y)@f$ and radius @f$\rho@f$. The
transverse momentum can be related to available quantities like


@f[
p_\mathrm{T} \propto \cdot q B \rho
@f]


with the charge @f$q@f$ and the magnetic field @f$B@f$. An intersection
between the straight line in the @f$rz@f$-plane with the @f$z@f$-axis gives an
estimate of the longitudinal impact parameter.
An illustration of seeds in the transverse plane is found in
the @ref seeding_figure "figure" below. Note that seeds can incorporate hits spread across all of
the layers shown, although this can be a configuration parameter.


@anchor seeding_figure

![Sketch of seeds in the transverse plane for a number of tracks on four layers. Seeds can combine hits on any three of these layers. The shown seeds appear compatible with having originated in the center of the detector, which is also drawn.](tracking/seeding.svg) {width=50%}

# Track finding and track fitting {#track-finding-and-track-fitting}

In the track seeding and following approach, track candidates are built from
the initial seeds. One method implemented in ACTS, the @ref
Acts::CombinatorialKalmanFilter "Combinatorial Kalman Filter" (CKF), uses the
*Kalman formalism*. Originally developed for monitoring and steering mechanical
systems, it can also be used to iteratively calculate a track estimate. After a
set of track candidates has been assembled and filtered (see @ref
ambiguity-resolution), an additional track fit is usually performed to extract
the best estimate of the track. The Kalman formalism can also be used for this,
with the addition of a smoothing step that has certain benefits.  Other fit
strategies exist, such as a global @f$\chi^2@f$ fit that minimizes the
distances between track-sensor intersections and measurements on all sensors at
the same time. One drawback of this method is the necessity to invert very
large matrices, which is computationally expensive.

In a track fit, the Kalman formalism can be shown to yield optimal estimates
for Gaussian uncertainties. This assumption breaks down when effects like
bremsstrahlung come into play. An extension of the Kalman Filter (KF) exists
that relies on the individual propagation of a set of trajectories, instead of
a single one, to model these biased uncertainties by a sum of Gaussian
components. This @ref Acts::GaussianSumFitter "Gaussian Sum Filter" (GSF) achieves better results when
fitting particles such as electrons, likely to undergo bremsstrahlung, and is
deployed in e.g. the ATLAS tracking chain.


## Kalman formalism and Kalman track fitter {#kalman-formalism}

> [!tip]
> See @ref Acts::KalmanFitter for documentation of the implementation of the
> Kalman Filter in the core library.

The basis of the Kalman formalism is a state vector, that can be identified
with the set of track parameters @f$\vec x@f$. Note that the concrete
parametrization plays a subordinate role in this context. Rather than building
an estimate of the state of a system in real time, a Kalman track fit can be
understood as estimating the parameters iteratively in steps. In the track
fitting application, each step is defined by a measurement to be included.
The evolution of the state vector is described by


@f[
  \vec x_k = \mathbf F_{k-1} \vec x_{k-1} + \vec w_{k-1},
@f]


where the linear function @f$\mathbf F_{k-1}@f$ transports the state vector at
step @f$k-1@f$ to step @f$k@f$. @f$\vec w_{k-1}@f$ is additional so-called process noise
that affects the transport additively. Each step has an associated
measurement, with the fixed relationship between the measurement and the state vector


@f[
  \vec m_k = \mathbf H_k \vec x_k + \epsilon_k.
@f]


Here, @f$\mathbf H_k@f$ is the *measurement mapping function*, which
transforms the state vector into the measurement space. In the ideal case,
this purpose can be achieved by a simple projection matrix, which extracts a
subspace of the state vector. Additionally, @f$\epsilon_k@f$ represents the
measurement uncertainty.

The Kalman fit process is divided into different phases:

1. **Prediction** of the state vector at the next step @f$k+1@f$ based on the information at the current step @f$k@f$.
2. **Filtering** of the prediction by incorporating the measurement associated to the step
3. **Smoothing** of the state vector by walking back the steps and using information for the subsequent step @f$k+1@f$ to improve the estimate at the current step @f$k@f$.

An illustration of these concepts is found in the @ref fig_kalman_filter "figure" below. Here,
a series of three sensors is shown with measurements on them. The KF
then predicts the track parameters at an intersection, shown in blue.
Subsequently, a filtered set of parameters is calculated as a mixture between
the measurement and the prediction. Not shown in this picture is the
smoothing step.

@anchor fig_kalman_filter

![Illustration of the KF. Two of the three stages, the prediction and the filtering are shown. The filtering updates the prediction with information from the measurement.](tracking/kalman.svg) {width=70%}


In many cases, the first two phases run in tandem, with prediction and
filtering happening alternatingly at each step. The smoothing phase is
launched once the last measurement has been encountered.
Starting from a state @f$k@f$, first, a prediction of the state vector at the
next measurement location is obtained via

@anchor kf_pred

@f[
  \vec x_k^{k-1} = \mathbf F_{k-1} \vec x_{k-1},
@f]


with the linear transport function from above. @f$\vec x_k^{k-1}@f$ is
the prediction of the state vector at step @f$k@f$ based on step @f$k-1@f$. The next
stage is the filtering. Here, the state vector is refined by taking into
account the measurement at the current step. Following one of two variants of
filtering from @cite Fruhwirth:1987fm, the gain matrix formalism, the state
vector is updated like


@f[
  \vec x_k = \vec x_k^{k-1} + \mathbf K_k \left( \vec m_k - \mathbf H_k \vec x_k^{k-1} \right),
@f]


with the *Kalman gain matrix*


@f[
  \mathbf K_k = \mathbf C_k^{k-1} \mathbf H_k^\mathrm{T}
    \left(
      \mathbf V_k + \mathbf H_k \mathbf C_k^{k-1} \mathbf H_k^\mathrm{T}
    \right)^{-1}
    .
@f]


Note that @f$\vec x_k@f$ is the filtered state vector at step @f$k@f$,
based on information from previous steps and step @f$k@f$ itself. This is in
contrast to @f$\vec x_k^{k-1}@f$, which is the prediction of the state vector at
step @f$k@f$ based on @f$k-1@f$, and is used to calculate the filtered state vector.
One input to these equations is the covariance matrix prediction @f$\mathbf
C_k^{k-1}@f$ at step @f$k@f$ based on step @f$k-1@f$, which can be written as


@f[
  \mathbf C_k^{k-1}  = \mathbf F_{k-1} \mathbf C_{k-1} \mathbf F_{k-1}^\mathrm{T} + \mathbf Q_{k-1}
@f]


in the linear version from @cite Fruhwirth:1987fm, with the
covariance @f$\mathbf C_{k-1}@f$ at step @f$k-1@f$, and the covariance @f$\mathbf
Q_{k-1}@f$ associated with @f$\vec w_{k-1}@f$ from above. Also needed is @f$\mathbf
V_k@f$, which is the covariance associated with @f$\epsilon_k@f$, effectively
representing the measurement uncertainty.

Similar to the state vector itself, the corresponding covariance matrix is
also filtered using

@anchor kf_cov_pred

@f[
  \mathbf C_k = \left( \mathbb 1 - \mathbf K_k \mathbf H_k \right) \mathbf C_k^{k-1}.
@f]


In the smoothing phase, the state vector at step @f$k@f$ is improved using the
information from the subsequent step @f$k+1@f$ using


@f[
  \vec x_k^n = \vec x_k + \mathbf A_k \left( \vec x_{k+1}^n - \vec x_{k+1}^k \right).
@f]


Here, @f$\vec x_{k+1}^n@f$ is the smoothed state vector and @f$\vec
x_{k+1}^k@f$ the predicted state vector at the subsequent step @f$k+1@f$. Also
needed is the *smoother gain matrix*


@f[
  \mathbf A_k = \mathbf C_k \mathbf F_k^\mathrm{T} \left( \mathbf C^k_{k+1} \right)^{-1},
@f]


with the predicted covariance at step @f$k+1@f$, @f$\mathbf C_k^{k+1}@f$.
Finally, the covariance at the current step @f$k@f$ can also be smoothed with


@f[
  \mathbf C_k^n = \mathbf C_k + \mathbf A_k \left(\mathbf C_{k+1}^n - \mathbf C_{k+1}^k \right) \mathbf A_k^\mathrm{T}.
@f]


After smoothing, the parameter estimate at the first step contains information
from all other measurement states. As mentioned above, in case the
uncertainties entering the Kalman fit are Gaussian distributions without
biases, the KF can be shown to be the optimal solution minimizing mean
square estimation error. However, certain caveats exist. The KF assumes
that a linear transport function @f$\mathbf F@f$ exists that can propagate the
state vector. In the presence of inhomogeneous magnetic fields this is not
the case. Instead of explicitly applying @f$\mathbf F@f$ to the state vector for
the prediction, the ACTS KF turns to the numerical integration,
discussed in @ref numerical-integration. With it, the prediction from
@ref kf_pred "this equation" is simply the intersection of the extrapolated trajectory
with the next sensitive surface. Aside from this, @f$\mathbf F@f$ is also used to
transport the covariance between steps (see @ref kf_cov_pred "here"). Here, the
semi-analytical method for covariance transport in the numerical integration
can be used. @f$\mathbf F@f$ can then be identified with the transport
Jacobian accumulated between surfaces.

For smoothing, two possibilities exist to obtain the needed covariances from
subsequent measurement steps. Either, the inverse transport Jacobian is used
and applied, in a way similar to @ref kf_cov_pred "this equation", or the numerical
integration is executed again in an inverse fashion, propagating from the
subsequent step to the current one.

## Combinatorial Kalman Filter {#combinatorial-kalman-filter}

> [!tip]
> See @ref Acts::CombinatorialKalmanFilter for information on the CKF
> implementation found in the core library.

As mentioned above, the Kalman formalism can be used for track finding. In this
case, the smoothing step can be dropped, as the resulting track candidates are
likely to be refit regardless, therefore saving some time. The CKF explores the
event starting from an initial track seed. It does this by considering not only
a single sequence of measurements, but allowing the branching of the fit at
each sensitive surface that is encountered. To this end, all or a subset of
measurements that are found on each surface are considered. Measurements are
selected based on their compatibility with the current state estimate, by using
their residuals. A predicted residual


@f[
  \vec r_k^{k-1} = \vec m_k - \mathbf H_k \vec x_k^{k-1},
@f]


and a filtered residual


@f[
  \vec r_k = \vec m_k - \mathbf H_k \vec x_k,
@f]


can be defined, depending on which state estimate is compared with
the measurement @f$\vec m_k@f$. Using the filtered residual, an effective
@f$\chi^2@f$ increment


@f[
  \chi^2_+ = \vec r_k^\mathrm{T}
  \left[ \left( \mathbb 1 - \mathbf H_k  \mathbf K_k \right)  \mathbf V_k \right]^{-1}
  \vec r_k
@f]


can be calculated. The global @f$\chi^2@f$ of the track candidate can
be calculated as the sum of all @f$\chi^2_+@f$ across the steps. Measurements
with a large @f$\chi^2_+@f$ are considered as outliers, which have low
compatibility with the trajectory. By branching out for measurements below a
certain @f$\chi^2_+@f$, and following the branches, a tree-like structure of
compatible track candidates originating from a track seed is assembled. This
feature is shown in the @ref fig_tracking_ckf "figure" below, which displays a circular
trajectory, and a set of iteratively assembled track candidates. Basic
quality criteria can be applied at this stage, to remove bad candidates. A
dedicated @ref ambiguity-resolution.
selects the candidates most likely to belong to real particle tracks.

@anchor fig_tracking_ckf

![Illustration of the way the CKF iteratively explores measurements from a seed outwards. Measurements are added successively, and can be shared between the resulting track candidates. Shown in green is a circular *real* trajectory.](tracking/finding.svg) {width=50%}

# Ambiguity resolution {#ambiguity-resolution}

Due to the combinatorial nature of track finding, and to achieve high
efficiencies, this set of candidates is often large, and contains a
non-negligible fraction of *fake* candidates. These fake candidates are either
completely combinatorial, or arise from real particle measurements with
combinatorial additions. Track candidates coming from a single seed necessarily
share a common stem of measurements. Measurements can potentially also be
shared between candidates from different seeds.

One possibility to resolve this (as is done in e.g. ATLAS tracking) is an
ambiguity resolution algorithm, that attempts to filter out as many undesirable
tracks as possible. This is implemented by means of a scoring function, that
combines properties of the track parameters. Higher scores are correlated with
a larger probability to be a desirable track candidate. A larger number of hits
results in an increase in the score, as longer compatible hit chains are less
likely to be random combinations. On the other hand, missing hits in sensors
where a hit was expected negatively impact the score.  Experiment specific
scoring of hits from different subsystems is also implemented. The overall
@f$\chi^2@f$ value computed for the track candidate also plays a role. Candidates
that share hits with other candidates are penalized. Another quantity is the
measured particle @f$p_\mathrm{Y}@f$, which enters the score, to give preference to
tracks with large momenta. For tracks containing measurements with a
substantial local @f$\chi^2_+@f$ at the start or end of the trajectory, the
ambiguity resolution step can also attempt to remove these hits, and determine
whether a refit without them yields a more favorable global @f$\chi^2@f$.

Finally, the output of the ambiguity resolution step is a set of track candidates
that contain an enhanced fraction of tracks from actual particles, while fake
tracks are suppressed. They are passed into the final precision fit outlined
in @ref track-finding-and-track-fitting, to extract the parameter estimate, and used
for further aspects of reconstruction.

# Vertex reconstruction {#vertex-reconstruction}

> [!tip]
> See @ref vertexing for more information on the vertexing as implemented in ACTS.

A vertex is a point within the detector, where an interaction or a
decay occurred. We distinguish between primary vertices (from
collisions/interactions) and secondary vertices (from subsequent particle
decays), see the @ref vertexing_illust "figure" below. Primary vertices are further divided
into hard-scatter and pile-up vertices. While primary vertices are located in
the luminous region, secondary vertices are slightly displaced due to the finite
 life time of the decaying particle.

@anchor vertexing_illust

![Illustration of a set of three vertices in a proton-proton collision. We distinguish between primary hard-scatter, primary pile-up, and secondary vertices.](tracking/vertexing.svg) {width=60%}

Vertices play an important role in higher-level reconstruction algorithms. For
example, secondary vertices can help with the identification of particles:
During *b-tagging*, a displaced vertex located inside a jet is a sign for the
decay of a @f$b@f$-hadron.

In analogy to track reconstruction, vertex reconstruction can be divided into
two stages: vertex finding and vertex fitting. As a first step of vertex
finding, we compute a rough estimate of the vertex position from a set of
tracks. This first estimate can be calculated in many different ways, and is
referred to as "vertex seed". Seeding algorithms differ for primary and
secondary vertexing. For primary vertex seeding, one option is to use a
histogram approach to cluster tracks on the @f$z@f$-axis @cite Piacquadio_2010.
This is based on the assumption that primary vertices will be close to the
beamline. Other approaches model tracks as multivariate Gaussian distributions
and identify regions of high track density as vertex seeds @cite schlag_2022.
For secondary vertexing, seeds are formed from pairs of reconstructed tracks as
the constraint to the beamline does not apply.

Once a vertex seed is determined, tracks that are compatible with it are
selected as part of the vertex finding.

Before the vertex fit, we linearize tracks in the vicinity of the vertex seed
under assuming that they follow a helical (for constant magnetic field) or
straight (for no magnetic field) trajectory @cite Piacquadio_2010. The vertex
fitter then uses this linearization to improve the position of the vertex seed.
Furthermore, the track momenta are refitted under the assumption that the tracks
originate at the vertex @cite Fruhwirth:1987fm @cite Billoir:1992yq .

One issue with an approach like this is that the assignment of tracks to
vertices is ambiguous. As an improvement, one can perform a multi-vertex fit,
where vertices compete for tracks. This means that one track can be assigned to
several vertices. Their contribution to each vertex fit is determined by a
weight factor, which, in turn, depends on the tracks' compatibility with respect
to all vertices @cite Fruhwirth:2005.

A flowchart of a multi-vertex reconstruction chain is shown in
the @ref vertexing_flowchart "figure" below.

@anchor vertexing_flowchart

![Simplified flowchart of multi-vertex reconstruction. From a set of seed tracks, we first compute a rough estimate of the vertex position, i.e., the vertex seed. Then, we evaluate the compatibility of all tracks with the latter. If a track is deemed compatible, it is assigned a weight and attached to the vertex seed. Next, the vertex seed and all previously found vertices that share tracks with it are (re-)fitted. Finally, after convergence of the fit, we check whether the vertex candidate is merged with other vertices and discard it if that is the case. For the next iteration, all tracks that were assigned to the vertex seed and that have a weight above a certain threshold are removed from the seed tracks.](tracking/vertexing_flowchart.svg)
