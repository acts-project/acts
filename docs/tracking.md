# Tracking in a nutshell

Track reconstruction is the process to recover the properties of a charged
particle from a set of measurements caused by interaction with some form of
sensitive detector. The goal is to find which measurements are likely to have
been caused by which particle, group them accordingly, and estimate the
associated trajectory. Such charged particle trajectories form the basic
input to the majority of higher-level reconstruction procedures in many
cases.

:::{figure} /figures/tracking/tracking.svg
:width: 400px
:align: center

Illustration of a track reconstruction chain starting from spacepoints to fully
formed tracks
:::

This section provides a high-level view of a track reconstruction chain, and is largely based on [^phd_paul].

## Charged particle detection

The first step in the chain to reconstruct charged particles is their
detection using sensitive elements. Charged particle detection can be
achieved in a variety of ways with very different technologies. It can be
achieved by measuring the interaction of charged particles with matter. When
this occurs, the interacting particle typically ionizes the surrounding
material. Particle detectors make use of this fact by converting the
resulting charge into a measurable signal in various ways.

(segmentation)=
:::{figure} /figures/tracking/segmentation.svg
:align: center
:width: 400px

Illustration of one-dimensional (a) and two-dimensional segmentation (b) of a silicon sensor.
:::

A very common electronic detection approach is the use of semiconducting
particle detectors, often made of silicon. When a charged particle traverses
such a sensor, it ionizes the material in the depletion zone, caused by the
interface of two different semiconducting materials. The result are pairs of
opposite charges. These charge pairs are separated by an electric field and
drift toward the electrodes. At this point, an electric signal is created
which can be amplified and read out. By means of segmentation, the measured
signal can be associated with a location on the sensor. Silicon sensors are
usually segmented in one dimension (*strips*) or in two dimensions
(*pixels*) (see {numref}`segmentation`).

(track_parametrization)=

## Track parametrization

To express the properties of a particle's trajectory, a choice of parameters
has to be made. The parameters need to express all the relevant quantities of
interest. In the presence of a magnetic field, which affects charged
trajectories, the global position and momentum, as well as the charge are
needed to fully specify the particle properties. In addition, a time parameter
can be included. Apart from the global reference frame, track quantities often
need to be represented with respect to a surface. This can be achieved with a
parametrization like

\begin{equation*}
  \vec x = \left(l_0, l_1, \phi, \theta, q/p, t\right)^T
\end{equation*}

although other parameter conventions exist as well.
{numref}`parameters` illustrates this choice of parameters. $l_0$, $l_1$ are
the local coordinates of the corresponding surface, $\phi \in [-\pi,\pi)$ and
$\theta \in [0,\pi]$ are the angles in the transverse and longitudinal
direction of the global frame, expressed with respect to the current location
along the trajectory, as indicated in {numref}`parameters` (b). $\theta$ is
the polar angle relative to the positive $z$-axis, and $\phi$ is the azimuth
angle in the transverse plane. Finally, $q/p$ combines the charge of the
particle with the inverse momentum. In {numref}`parameters` (a), the global
momentum vector $\vec p$ is shown, which can be recovered from the parameters
$\vec x$ using $\phi$, $\theta$ and $q/p$.

(parameters)=
:::{figure} /figures/tracking/parameters.svg
:width: 500px
:align: center

Illustration of the parametrization of a particle track with respect to a
two-dimensional surface. (a) shows the local position, global momentum and
their corresponding uncertainties. (b) displays the angles $\phi$ and $\theta$
in the transverse and longitudinal planes.

:::

(perigee)=
:::{figure} /figures/tracking/perigee.svg
:width: 200px
:align: center

Illustration of the perigee parametrization which uses the point of closest
approach relative to a reference point. The impact parameter $d_0$, the
position $l$ and the momentum vector $\vec p$ are shown.

:::

Aside from the nominal quantities captured in $\vec x$, the related
uncertainties and correlations need to be taken into account as well. They
can be expressed as a $5\times 5$ covariance matrix like

\begin{equation*}
  C =
  \begin{bmatrix}
   \sigma^2(l_0)& \text{cov}(l_0,l_1) & \text{cov}(l_0, \phi) & \text{cov}(l_0, \theta) & \text{cov}(l_0, q/p) \\
   . & \sigma^2(l_1) & \text{cov}(l_1, \phi) & \text{cov}(l_1, \theta) & \text{cov}(l_1, q/p) \\
   . & . &  \sigma^2(\phi) & \text{cov}(\phi,\theta) & \text{cov}(\phi, q/p) \\
   . & . & . & \sigma^2(\theta) & \text{cov}(\theta, q/p) \\
   . & . & . & . & \sigma^2(q/p)
  \end{bmatrix}
\end{equation*}

Here, $\text{cov}(X,Y)$ is the covariance of variables $X$ and $Y$, while
$\sigma^2(X)$ are the regular variances. As the covariance matrix $C$ is
symmetric, only the upper right half is shown in the matrix above. The
uncertainties associated with the local position, as well as the momentum
direction are indicated in {numref}`parameters` (a) as an ellipse and a cone
around the momentum vector $\vec p$, respectively. 

(particle_propagation)=

## Particle propagation

A central for track reconstruction is the ability to calculate the trajectory
of a charged particle, given its properties at a given point. This process,
called *particle propagation* or *extrapolation*, is used to predict a
particle's properties after it has travelled a certain distance. In many cases,
the projected intersection with various types of surfaces is desired. The
trajectory of a charged particle is governed by the [magnetic
field](core/magnetic_field.rst) through which it travels, as well as any
[material effects](core/material.rst). In case of a homogeneous magnetic field,
and in the absence of material interaction, the particle follows a helical
trajectory. Such a helix can be calculated purely analytically, although
intersections require numerical methods nevertheless.

Often, magnetic fields are not homogeneous, however. In the presence of such
changing fields, the corresponding differential equations of motions need to be
solved using numerical integration techniques.

### Numerical integration

In ACTS, numerical integration is done using the *Runge-Kutta-Nyström* (RKN) method.
Commonly used in the variant at fourth order, it describes how to calculate a
solution to an initial value problem that can be formulated generically like

$$
  \frac{dy}{dt} = f(t,y), \qquad y(t_0) = y_0,
$$

 where $y_0$ refers to the initial value of $y$ at $t_0$, and
$f(t,y)$ is the functional form describing the dynamics. The method then
successively approximates the analytical solution $y(t)$ in a stepwise
fashion. At each step $(t_n, y_n)$, the goal is effectively to approximate
the next step $y(t_{n+1})$. Using a step size $h$, the algorithm evaluates
the function $f$ at four points $k_{1-4}$:

$$
  \begin{aligned}
    k_1 &= f(t_n, y_n) \\
    k_2 &= f\left( t_n + \frac h 2, y_n + h \frac{k_1} 2 \right) \\
    k_3 &= f\left( t_n + \frac h 2, y_n + h \frac{k_2} 2 \right)\\
    k_4 &= f\left( t_n + h, y_n + hk_3 \right).
  \end{aligned}
$$


(rk)=
:::{figure} /figures/tracking/rk.svg
:width: 400px
:align: center

Illustration of the RKN method approximating a first order
differential equation. Shown is the calculation of an estimate $y_{n+1}$ at
$t_{n+1} = t_n + h$, based on the current step $(t_n,y_n)$. Shown are the four
distinct points at which function $y(t)$ is evaluated, and which are blended to
form the estimate.

:::


Looking at {numref}`rk`, the meaning of these four points in relation
to the step size $h$ can be understood. $k_1$ is the derivative at the
current location, $k_{2,3}$ use $k_1$ and $k_2$ respectively to calculate two
envelope derivatives at $h/2$ and $k_4$ uses $k_3$ to make an estimate of the
derivative at $h$. Combining $k_{1-4}$, $(t_{n+1},y_{n+1})$ can be calculated
as the approximation of $y(t_{n+1})$ like

$$
  \begin{aligned}
    y_{n+1} &= y_n + \frac 1 6 h ( k_1 + 2 k_2 + 2 k_2 + k_4)\\
    t_{n+1} &= t_n + h
  \end{aligned}
$$

by effectively averaging the four derivatives. It is apparent that
the step size crucially influences the accuracy of the approximation. A large
step size weakens the approximation, especially if the magnetic field changes
strongly. On the other hand, a too small step size will negatively affect the
execution time of the algorithm.


The Runge-Kutta-Nyström method from above can be adapted to handle second order
differential equations, as is needed for the equations of motion in question,

$$
  \frac{d^2 \vec r}{ds^2} = \frac q p \left( \frac{d\vec r}{ds} \times \vec B (\vec r) \right) = f(s, \vec r, \vec T), \qquad \vec T \equiv \frac{d \vec r}{ds},
$$

with the global position $\vec r$, the path element $s$, the
normalized tangent vector $\vec T$ and the magnetic field $\vec B(\vec r)$ at
the global position. A slight modification of $k_{1-4}$ is also required,
incorporating the first derivative of $f(s, \vec r, \vec r')$, finally
leading to

$$
  \begin{aligned}
    \vec T_{n+1} &= \vec T_n + \frac h 6 (k_1 + 2k_2 + 2k_3 + k_4) \\
    \vec r_{n+1} &= \vec r_n + h \vec T_n + \frac{h^2}{6} (k_1 + k_2 + k_3).
  \end{aligned}
$$

A strategy exists to dynamically adapt the step size according to the magnetic
field strength, with the definition of a target accuracy that the algorithm
tries to achieve. Here, the step size $h$ will successively be decreased and
the approximation recalculated until the accuracy goal is achieved. Even with
these additional calculations, the approach is still preferable over a
consistently low step size.

### Covariance transport

Aside from the prediction of the track parameters at a given path length, a key
ingredient to many dependent applications are the uncertainties in the form of
the associated covariance matrix $C$. Conversions between covariance matrices
$C^i\to C^f$ can generally be achieved like

$$
  C^f = J \cdot C^i \cdot J^T,
$$

using the Jacobian matrix

$$
  J = \begin{bmatrix}
    \frac{\partial l_0^f}{\partial l_0^i} & \cdots & \frac{\partial l_0^f}{\partial (q/p)^i} \\
    \vdots & \ddots & \vdots \\
    \frac{\partial (q/p)^f}{\partial l_0^i} & \cdots & \frac{\partial (q/p)^f}{\partial (q/p)^i}
  \end{bmatrix},
$$

between initial and final parameters $\vec x^i$ and $\vec x^f$. The
task therefore becomes calculating the necessary Jacobians to achieve correct
transformation.

One part is the transformation between different
coordinate systems, but at the same location along the trajectory. For this
purpose, generic Jacobians can be calculated between each coordinate system
type, and a common coordinate system. The common coordinate system used for
this purpose is the curvilinear frame, which consists of the global direction
angles, and a plane surface located at the track location, with its normal
aligned with the track momentum. By using Jacobians to the curvilinear frame
and the corresponding inverse matrices, conversions between any two
coordinate systems can be performed.

The second part is the calculation of the evolution of the covariance matrix
during the propagation between surfaces. To this end, a semi-analytical
method which calculates the effective derivatives between two consecutive
RKN steps can be used. By accumulating the Jacobian
matrices calculated for each step, the effective Jacobian between the
starting point and the destination can be obtained.

## Material effects

Charged particles interact with matter as they pass through it. Since
particle detectors inevitably consist of some form or material, this effect
cannot be completely avoided. By building tracking detectors as light as
possible, and arranging passive components, such as services and support
structures carefully, the material a particle encounters before being
measured can be reduced. Charged particles traversing any form of matter
undergo elastic and inelastic interactions with the atomic structure of the
material, depending on the particle properties.

(multiple_scattering)=

:::{figure} /figures/tracking/multiple_scattering.svg
:width: 200px
:align: center

Illustration of the effect of multiple scattering on the
trajectory of a charged particle passing through a block of material.
Entering from the left, it undergoes a series of scattering events,
deflecting the trajectory statistically, before exiting on the right

:::

In elastic interactions, the particle does not lose a significant amount of
energy, while its trajectory is affected. {numref}`multiple_scattering` shows a
sketch of the way multiple Coulomb scattering affects the direction of a
particle trajectory. In addition, a shift in the transverse plane relative to
the incident direction can occur. As the scattering events occur in
statistically independent directions, the means of both the deflection and
offset tends toward zero as the number of scatters increases. Therefore, in the
numerical particle propagation, this can be accounted for by simply increasing
the uncertainties associated with the direction, depending on the amount of
material encountered.

On the other hand, there are interactions during which the particle loses
some of its energy. Relevant processes here are ionization, as well as
bremsstrahlung for light particles like electrons. For hadronic particles,
hadronic interactions with the nuclei of surrounding material is another
process of interest. In such hadronic interactions, the incoming particle
often disintegrates, and does not propagate further. Since the size of
ionization losses only fluctuates to a small degree for thin layers of
material, they can usually be accounted for by reducing the trajectory energy
correspondingly. For bremsstrahlung, where fluctuations are much larger,
dedicated techniques are needed (see [](gsf)).

Two main approaches are implemented in ACTS. The first approximates the
material interaction by using a description that averages the real material
onto thin surfaces across the detector (more on this in
[](#geometry-and-material-modelling)). When the propagation encounters such a
surface, it retrieves the material properties, and executes parametrized
modifications to the particle properties and uncertainties. In the second
approach, material effects are continuously incorporated during propagation,
rather than at discrete locations. The latter approach is especially suited for
propagation through volumes of dense material, where the discretization of the
material distribution will not work as well.

## Geometry and material modelling

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
elements. {numref}`geometry_detail` (a) shows a sketch of such a detailed
geometry description. Shown as an example is a *layer* of silicon
sensors in a barrel configuration. The green rectangles represent the actual
sensitive surfaces, while other elements include cooling, readout and other components.

(geometry_detail)=

:::{figure} /figures/tracking/geometry_detail.svg
:align: center

Sketch of the way a fully detailed simulation geometry (a) models passive
elements, in addition to the sensitive elements shown in green. (b) shows a
simplified version, where all non-sensitive elements are approximated.

:::

In the majority of cases in track reconstruction, this detailed
geometry is unnecessary. During track reconstruction, the aforementioned
associated information needs to be accessible for measurements, so all
sensitive elements need to be included in some form. Passive elements, on the
other hand, are only required to factor in material interaction effects (see
[](#particle-propagation)). Moreover, the fully detailed geometry comes
at the disadvantage of introducing significant overhead during navigation. In
this process, an algorithm attempts to figure out which elements the particle
propagation needs to target, as the trajectory is likely to intersect them.
With a geometry description this precise, the navigation process becomes a
significant performance bottleneck.

(layer_barrel)=
:::{figure} /figures/tracking/layer_barrel.svg
:width: 300px
:align: center

Sketch of the way sensitive elements are grouped into layers.
Shown is an $xy$-view of a number of sensors, arranged as in e.g. the ATLAS
silicon detector barrels. The grouping is based on their mounting radius.
The layers are indicated in different colors.

:::

As a compromise between modelling accuracy and performance, ACTS
uses a simplified geometry model . It
focusses on the sensitive elements, which are strictly needed, while passive
elements are discarded from the explicit description and approximated.
{numref}`geometry_detail` (b) shows such a simplified geometry. Here, the
sensitive elements are still shown in green, and other elements are greyed
out, indicating that they are discarded. The sensitive elements are then
grouped into layers, as sketched in {numref}`layer_barrel`. How exactly the
grouping occurs depends on the concrete experiment geometrym. In some cases, the layers have the shape
of cylinder surfaces with increasing radii. This example is shown in the
figure in the transverse plane at radii $r_{1,2,3}$. In the endcaps, where
modules are arranged on disks, these are used as the layer shape. An
illustration of endcap disk layers can be found in {numref}`layer_ec`, where
six disks are located at six distinct positions in $\pm z_{1,2,3}$, and
shown in different colors.

(layer_ec)=
:::{figure} /figures/tracking/layer_ec.svg
:width: 400px
:align: center

Sketch of the way sensitive elements are grouped into layers.
Shown is a view of a number of sensors, arranged as in e.g. the ATLAS
silicon detector endcaps. They are grouped into disks based on their
mounting position in $z$. The layers are indicated in different colors.

:::

During particle propagation, the navigation makes use of this layer
system. Each layer contains a binned structure, which maps a bin to a set
of sensitive surfaces that overlap with the bin area. This is illustrated in
\cref{fig:geo_binning}, where the left picture shows the sensitive surface
structure of an exemplary endcap disk. The picture on the right overlays the
binning structure that can be used to enable fast retrieval of compatible
sensitive surfaces. By performing a simple bin lookup, the navigation can
ascertain which sensors it needs to attempt propagation to.


:::{figure} /figures/tracking/surface_array.svg
:width: 400px
:align: center

Illustration of the binning structure that is used to subdivide
layer surfaces. (a) shows two sensor rings of
different radii grouped into one disk layer. (b)
overlays the binning structure that the navigation queries for compatible
surfaces.

:::

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
\emph{approach surfaces} on either side. Their distance can be interpreted as
the thickness of the layer in question. Examples of these approach surfaces
can be found in {numref}`geometry_detail`, at the inner and outer radius.
Approach surfaces, and the boundary surfaces between volumes mentioned before,
are candidates to receive a projection of the surrounding material.
Additional artificial material layers can also be inserted to receive
projected material.

The projection procedure works by extrapolating test particles using the
fully detailed simulation geometry. During the extrapolation, the material
properties of the geometry are sampled in small intervals. Subsequently, the
same test particle is extrapolated through the tracking geometry. All
material samples are then assigned and projected onto the closest material
surface. Finally, the projection is averaged. The exact number and placement
of the material surfaces has to be optimized to yield a sufficiently accurate
representation of the inactive material in the detector.

The numerical integration uses these projected material surfaces. Whenever
such a surface is encountered in the propagation, the material properties are
retrieved, and the corresponding modifications to the trajectory are
executed. In case material is supposed to be integrated in a continuous way
(as mentioned in [](#particle-propagation)), volumes can also store an
effective volumetric material composition, which is queried by the numerical
integration when needed. As the actual physical location of the detection
hardware can vary over time, possible misalignment of the sensors needs to be
handled correctly (see [](#alignment-of-particle-sensors)).

## Clusterization

## Spacepoint formation and seeding

## Track finding and track fitting

## Ambiguity resolution

## Vertex reconstruction

## Alignment of particle sensors


[^phd_paul]: Gessinger-Befurt, Paul, 30.04.2021. Development and improvement of track reconstruction software and search for disappearing tracks with the ATLAS experiment. [10.25358/openscience-5901](https://doi.org/10.25358/openscience-5901)
