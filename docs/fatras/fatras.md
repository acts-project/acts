# Fatras

Fatras, short for `Fast Track Simulation` is monte carlo based track simulation that used the ACTS modules together with Monte Carlo driven integration of material effects. The original concept of Fatras stems from
the ATLAS experiment and is described [here](https://cds.cern.ch/record/1091969?ln=de). Fatras has been used
to simulate the data set of the [Tracking machine leearning challenge](https://arxiv.org/abs/1904.06778).

Fatras consists of:
 * an interface to particle generation/input module
 * an simulation kernel module
 * an geometric and smearing digitization module

## Principle

Fatras makes use of the geometry, navigation and propagation modules of ACTS. Given that the navigation module
of ACTS is capable of predicting the particle trajectory through the detector geometry, it can be used - together with the transport/stepper code of the propagator module as a particle simulator through the detector.

For Fatras, however, the material interaction of the particle is then changed to run using random numbers that are used for emulating physics processes. Technically, this is enabled with the `actor` mechanism of the propagator.


```{note}
Being a fast simulation, Fatras does not target to be as accurate as a full simulation setup, e.g. using the Geant4 toolkit. Its main purpose is to provide a fast, yet accurate enough simulation engine to support efficient algorithm development and large scale reconstruction benchmarking.
```

## Particle Input / Simulation

The Fatras module itself does not ship with a particle gun nor an interface to an event generator, but requires
that generated particles are translated into an internal `Particle` object.

:::{doxygenfunction} ActsFatras::Particle::Particle(Barcode particleId, Acts::PdgParticle pdg, Scalar charge, Scalar mass)
:::

```{note}
In the `Examples/Fatras` one can find an example how to fill the particle containers from e.g. a particle gun or
from Pythia8 input. Furthermore, a translation from `HepMC3` is shown.
```


## Particle Simulation

Charged and neutral particles are handled by the Fatras simulation kernel, via two dedicated simulators. Dedicated particle selectors can be configured to select which particles are then further processed within the `Simulation` kernel.

Once a particle is accepted for further processing, it is propagated with an appropriated `Propagator` instance 
through the detector geometry. Whenever detector material is passed, physics processes are evoked respecting particle type and kinematics, but also a given configuration. When a charged particle traverses sensitive detector material (i.e. in Acts terms a surface with sensitive definition)


```{warning}
Fatras uses the strict surface based model of the Acts TrackingGeometry and hence hits created with sensitive detectors are intersections with surfaces only, and hence do not per se reflect the concept of an energy deposit
within a certain path length, as e.g. known from Geant4. Such an expansion into a three dimensional hit, however, can be done with the (geometric) digitizaton module as described below.
```

### Physics Processes

The following physics processes are available and implement:

| Physics Process | Implementation | Status |
|-------|--------------------------|---------|
| Decay | The particle decay is outsourced to Geant4 | Off |
| Ionization Loss | Described by the Bethe-Bloch formula | On |
| Radiation Loss | Gescribed by the Bethe-Heitler formula | On |
| Photon conversion | Ad-hoc implementation | On |
| Multiple Scattering | Implemented as Gaussian mixture, general mixture, or Highland scattering | On |
| Nuclear interaction | A parameterized version exists | Off |

```{note}
Fatras is still under development and not all physics processes are fully validated, some are switched off by default. Modules will become available and switched on when validated against Geant4.
```

## Digitization

Fatras ships with a digitization module that allows to emulate both the readout and intrinsic uncertainty of a detector. Currently, there are two main digitization modes available:
 * parametric smearing, short *smearing* digitization
 * geometric channelization, also referred to as *geometric* digitization


 ### Smearing digitization

 This allows for the fastest and most simple approximation of detector readout, the smearing for the truth hit information with a given smearing function. Hereby it makes use of the fact that the measurement frame is assumed to be defined by the surface coordinate system and hence the bound track parameter representation provided by the propagator can be directly used for creating output measurements.

 The smear function has to follow the following API:

 ```c++
/// Smearing function definition for single track parameters.
///
/// The function takes the unsmeared parameter and returns the smeared value and
/// a standard deviation.
///
/// @tparam generator_t The type of the random generator.
template <typename generator_t>
using SingleParameterSmearFunction =
    std::function<Acts::Result<std::pair<double, double>>(double,
                                                          generator_t&)>;
 ```

And is then used in the bound parameter smearer:

:::{doxygenstruct} ActsFatras::BoundParametersSmearer
:::

```{note}
The smearing implementation assumes uncorrelated measurements between the variable sized bound parameters.
```
 ### Geometric digitization

For planar detector types, there exists a more realistic digitization that emulates the readout structure, i.e. in case of a pixel detector a pixel segmentation is used to evaluate how many (and which) pixels were traversed by the particle trajectory.

