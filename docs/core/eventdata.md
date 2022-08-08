# Event data

:::{attention}
This section is largely **outdated** and will be replaced in the future.
:::

## Track parametrization

:::{tip}
An introduction to track parametrization concepts can be found in
[](track_parametrization).
:::

A trajectory in a magnetic field is generally parameterized by a set of at least
five parameters (when being bound to a surface). Two different categories are
used in Acts: so-called bound parameters, i.e. parameter bound to a surface and
curvilinear parameters. Curvilinear parameters are defined in an implicit planar
and normal (to the track) reference plane that follows the track. The center of
the implicitly defined plane is at the current track position, while the normal
vector points along the momentum direction. Per definition the local parameters
of a curvilinear parametrization are thus fixed to (0,0).

In Acts the parametrization is can be changed, provided the according
transformations into global coordinates are given, it can be changed by adapting
the relevant `coordinate_transformation` definition. This shows an excerpt of
the default implementation:

```cpp
namespace Acts {
namespace detail {
  /**
   * @brief helper structure summarizing coordinate transformations
   */
  struct coordinate_transformation
  {
    typedef ActsVector<ActsScalar, Acts::NGlobalPars> ParVector_t;

    static Vector3
    parameters2globalPosition(const ParVector_t& pars, const Surface& s)
    {
      return s.localToGlobal(Vector2(pars(Acts::eBoundLoc0), pars(Acts::eBoundLoc1)),
                             parameters2globalMomentum(pars));
    }

    static Vector3
    parameters2globalMomentum(const ParVector_t& pars)
    {
      Vector3 momentum;
      double         p     = std::abs(1. / pars(Acts::eBoundQOverP));
      double         phi   = pars(Acts::eBoundPhi);
      double         theta = pars(Acts::eBoundTheta);
      momentum << p * sin(theta) * cos(phi), p * sin(theta) * sin(phi),
          p * cos(theta);

      return momentum;
    }

    static ParVector_t
    global2curvilinear(const Vector3&,
                       const Vector3& mom,
                       double                charge)
    {
      ParVector_t parameters;
      parameters << 0, 0, mom.phi(), mom.theta(),
          ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.mag();

      return parameters;
    }
  // ...
  };
} // end of namespace detail
} // end of namespace Acts
```

Changing the default parameter definition and transformation needs recompilation
of the Acts Software and redefinition of the relevant plugin:
 
```cmake
target_compile_definitions (Acts::Core PUBLIC -DACTS_PARAMETER_DEFINITIONS_PLUGIN="${ACTS_PARAMETER_DEFINITIONS_PLUGIN}")
```

### Neutral and charged representations

Track parameters come in a charged and neutral flavor, the flavor is defined using a template parameter,
which either is a `ChargedPolicy` class for charged track parameter representation:

```cpp
namespace Acts {
    typedef SingleTrackParameters<ChargedPolicy>            TrackParameters;
    typedef SingleCurvilinearTrackParameters<ChargedPolicy> CurvilinearTrackParameters;
    typedef SingleBoundTrackParameters<ChargedPolicy>       BoundTrackParameters;
}  // end of namespace Acts
```

Or, respectively, a `NeutralPolicy` object

```cpp
namespace Acts {
    typedef SingleTrackParameters<NeutralPolicy> NeutralParameters;
    typedef SingleCurvilinearTrackParameters<NeutralPolicy>
                                              NeutralCurvilinearTrackParameters;
    typedef SingleBoundTrackParameters<NeutralPolicy> NeutralBoundTrackParameters;
} // end of namespace Acts  
```

### Multi-variant representation

Multi-variant fitters, such as the Gaussian Sum filter rely on a multi-component
description of the track, which requires the definition of multi-variant track
parameters. For Acts, as Propagator and Extrapolator are written for type
templates, a multi-variant definition of track parametrization is *planned* in
order to integrate directly with the core software.

## Measurements

Measurements exist as *uncalibrated* and *calibrated* types. While
*uncalibrated* measurements are the direct output of the data formation, during
the track fit, when trajectory information (or information about the trajectory
hypothesis) is available, certain calibration steps can help to improve the
track fit quality.

### Uncalibrated measurements

Measurements in a detector can be one to many-dimensional (covering the full
range up to the full track parametrization).

```cpp
template <typename Identifier, BoundIndices... params>
class Measurement
{
  // check type conditions
  static_assert(std::is_copy_constructible<Identifier>::value,
                "'Identifier' must be copy-constructible");
  static_assert(std::is_move_constructible<Identifier>::value,
                "'Identifier' must be move-constructible");
  static_assert(std::is_copy_assignable<Identifier>::value,
                "'Identifier' must be copy-assignable");
  static_assert(std::is_move_assignable<Identifier>::value,
                "'Identifier' must be move-assignable");

private:
  // private typedef's
  typedef ParameterSet<params...>
      ParSet_t;  ///< type of the underlying ParameterSet object
// ...
};
```
          
In order to minimize the computational cost (and differently from the original
ATLAS code base), the dimension of the measurement has to be fixed at compile
time.          

### Calibrated measurements

Calibrated measurements are temporary objects needed for track fitting and
belonging to a track. *Still to be defined*
