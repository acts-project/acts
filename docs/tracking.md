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

### Numeric integration

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

## Geometry and material modelling

## Clusterization

## Spacepoint formation and seeding

## Track finding and track fitting

## Ambiguity resolution

## Vertex reconstruction

## Alignment of particle sensors


[^phd_paul]: Gessinger-Befurt, Paul, 30.04.2021. Development and improvement of track reconstruction software and search for disappearing tracks with the ATLAS experiment. [10.25358/openscience-5901](https://doi.org/10.25358/openscience-5901)
