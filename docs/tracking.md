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
usually segmented in one dimension (**strips**) or in two dimensions
(**pixels**) (see {numref}`segmentation`).


## Track parametrization

In order to be able to express the properties of a particle's trajectory, a choice
of parameters has to be made. The parameters need to be able to express all
the relevant quantities of interest. In the presence of a magnetic field,
which affects charged trajectories, the global position and momentum, as well
as the charge are needed to fully specify the particle properties. In
addition, a time parameter can be included, but is not used in this
chapter. Apart from the global reference frame, track quantities often need
to be represented with respect to a surface. This can be achieved with a
parametrization like

\begin{equation*}
  \vec x = \left(l_0, l_1, \phi, \theta, q/p\right)^T
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

## Particle propagation
### Numeric integration
### Covariance transport

## Geometry and material modelling

## Clusterization

## Spacepoint formation and seeding

## Track finding and track fitting

## Ambiguity resolution

## Vertex reconstruction

## Alignment of particle sensors


[^phd_paul]: Gessinger-Befurt, Paul, 30.04.2021. Development and improvement of track reconstruction software and search for disappearing tracks with the ATLAS experiment. [10.25358/openscience-5901](https://doi.org/10.25358/openscience-5901)
