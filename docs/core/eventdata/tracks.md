:::{attention}
This section is **incomplete!**
:::

:::{contents}
:::

(edm_tracks)=
# High-level Track Event Data Model

Track information in ACTS can be divided into two parts: track-level
information and track state-level information.

Track-level information are properties that relate to the full track. This
includes the fitted track parameters with respect to some reference point,
often the origin of the detector or the beamspot. It can also include summary
information from the track finding stage, like the overall number of clusters
that were used in the creation of the track, or the fit quality from the track
fit.

Tracks are built-up from track states, where each track state corresponds to a
discrete state determining the track properties. This mainly includes
measurements states, expected intersections with sensors where no measurement
was found ({term}`holes<Hole>`), and intersections with known passive material.  The
{term}`EDM` allows building up a track from these track states iteratively. For
example, the Kalman Filter will append track states to the sequence whenever it
encounters a sensitive detector layer. The content of the track states is
defined such that the fitter can store all relevant information, with as little
need for extra information as possible. It is also designed to be flexible
enough to support different fitters, which might require different information
to be stored, as well as the combinatorial Kalman Filter, which produces a tree
of track states, instead of a fully linear sequence.

Ultimately, each output track is associated with a well-defined sequence of
track states, allowing downstream consumers of the {term}`EDM` to access the fully
detailed information produced during track reconstruction.

---

This section will first discuss the architecture and general implementation of
the Track {term}`EDM` in [](#edm_track_conceptual), and then report on the API
of the involved classes in [](#edm_track_api).

(edm_track_conceptual)=
## Conceptual

### Architecture

The Track {term}`EDM` is structured such that the memory-layout can be
{term}`SoA`, while presenting an object-oriented interface for convenient
usage.

{numref}`track_proxy` shows this object-oriented access model for
the example of the track container and track proxy object.  The track container
holds vectors of the various pieces of information, and has methods to add a
track, and to allow iteration over all tracks. This iteration, or index based
access, yields a track proxy object, which exposes the properties as methods
returning references, while internally only holding a pointer to and an index
into the track container. The types are built in a way that preserves
const-correctness, i.e. even though a track proxy is a value type which can be
copied, it will not allow modification of the underlying track container if it
is immutable:

```cpp
auto mutableTrackContainer = /*...*/;
auto trackProxy = mutableTrackContainer.getTrack(5); // is mutable
const auto& constTrackProxy = trackProxy; // is const
// ...
auto constTrackContainer = /*...*/;
auto trackProxy = trackContainer.getTrack(5); // is const, even as an lvalue
```

(track_proxy)=
:::{figure} figures/proxy.svg
:align: center
Illustration of the proxy pattern used in the track {term}`EDM`. The track
proxy logically represents a single track, and points to the data stored in
the track container.
:::

The track {term}`EDM` is fully agnostic to the concrete persistency framework
of an experiment. This avoids having to convert the data between different
representations, if implemented correctly.

### Implementation

To make the {term}`EDM` implementation independent of an experiment persistency
framework, it is separated into a *frontend layer* and a *backend layer*. The
frontend layer contains user-facing getters and setters, as well as convenience
methods that can be helpful. These methods are located either in the proxy
objects or in the containers, depending on whether they operate on a single
element or the entire container.

Overall, there are four main classes that make up the frontend layer:
{class}`Acts::TrackProxy`, `TrackContainer`, `TrackStateProxy` and `MultiTrajectory`.  The
latter serves as the track state container, where the name indicates that it is
able to handle a branching tree structure of track states.  `TrackProxy` and
`TrackStateProxy` expose methods to get the local track parameters and
covariance, corresponding reference surface, and also includes global
statistics like the total number of measurements, {term}`outliers<Outlier>` or
{term}`holes<Hole>` in case of `TrackProxy`.  `TrackProxy` also has a method to
conveniently iterate over the associated track states from the outside inwards,
yielding `TrackStateProxy` objects from the track state container.

In case of `TrackStateProxy`, functionality is exposed in the frontend
layer to allocate optional components, with the goal of reduced memory
footprint. There are two main uses of this: track parameters and measurements.
The track-state {term}`EDM` supports storing up to three sets of local track parameters
and covariance matrices, modeled after the information the Kalman Filter
formalism needs to store:

1. predicted parameter vector and covariance matrix
2. filtered parameter vector and covariance matrix
3. smoothed parameter vector and covariance matrix

In case of combinatorial track finding (see [](#track_finding), specifically
{numref}`tracking_ckf`), track hypothesis can start out with a common sequence
of track states, and then branch out when multiple compatible measurements are
encountered, as seen in {numref}`ckf_tree`.

The track state {term}`EDM` allows allocating only the track parameters that
are needed, and also allows sharing the same track parameters between multiple
track states, so that branching track states can share for example the same
predicted parameters. How this is achieved exactly is left to the backend
layer. Measurements are handled in a similar way, where the track finding
decides how much storage is needed based on the number of dimensions of an
incoming measurement. It then instructs the {term}`EDM` through the frontend
layer to ensure enough memory is available, where the specifics are again left
up to the backend layer.

The backend layer exposes an interface that is used by the frontend layer to
store and retrieve information.  It uses dedicated methods where needed, such as
for storing reference surfaces or source-link objects, which are lightweight
container objects for experiment-specific measurements. For the majority of
components, the frontend communicates with the backend through a single method
to obtain references to the underlying data. Components are accessed via hashes
of the component name, where the hashes are calculated at compile-time wherever
possible. The backend can then use the hashed component name to retrieve the
relevant memory. To allow directly manipulating the backing memory, the frontend
expects the backend to return references into the backing storage.

`TrackProxy` provides a method to copy a track between different track
containers, and only uses the frontend layer to accomplish this. This means that
copying tracks between different backend implementations is trivial.

(track_architecture)=
:::{figure} figures/edm_diagram.svg
:align: center
Diagram of the {term}`EDM` architecture. The frontend layer is used by other
ACTS components, and downstream clients. It is separated from the backend
layer by an interface. Conversion to and from EDM4hep is possible. Examples
of direct backend implementations are shown.
:::

{numref}`track_architecture` shows a diagram of the {term}`EDM` architecture. At the center
are the `TrackProxy` and `TrackContainer`. These classes are
produced by the track finding and track fitting components, and are the main
interface point with the clients of tracking.  In ACTS itself, all of the
performance monitoring and downstream reconstruction is either directly built on
top of these objects, or converts them into an internal {term}`EDM` on the use
case.  Behind the backend interface, the track container coordinates with both a
track state and a track backend, where a few examples are shown, and will be
discussed below.

(edm_track_api)=
## Track API

:::{doxygenclass} Acts::TrackProxy
:members:
:::

:::{doxygenclass} Acts::TrackContainer
:members:
:::

(edm_track_iteration)=
## Track state iteration and forward linking

By default, track states are connected as a one-directional linked list, where
each track state knows its *previous* track state. {numref}`ckf_tree` shows an
example of a track state tree, like it is constructed by the combinatorial
track finding.

In {numref}`ckf_tree` states $S_7$ and $S_6$ are the two {term}`tip states<tip
state>` of the track state tree, whil $S_1$ is the single {term}`stem state`.
In the case of combinatorial track finding starting from e.g. a
{term}`seed<Seed>`, it could be the location of the innermost {term}`space
point<Space point>`.

(ckf_tree)=
:::{figure} figures/ckf_tree.svg
:width: 400px
:align: center
Illustration of a branching multi-trajectory that is created during
combinatorial track finding.
:::

Each track object points at a single {term}`tip state` to define its track state sequence.
The {class}`Acts::TrackProxy` class has various methods to access the track state sequence:

```cpp
auto track = getTrackFromSomewhere();
for(const auto trackState : track.trackStatesReversed()) {
  // do something with track state
}

```

Note that {func}`Acts::TrackProxy::trackStatesReversed` iterates from the {term}`tip state` to
the {term}`stem state`, i.e. from the outside in.

:::{attention}
By-default, it is not possible to iterate *forward* through the track states on a track!
The track's track states need to be *forward-linked* for this to be possible.
:::

The reason for this is this:
As the trajectory branches at the second sensor into $S_2$/$S_3$, it is not
possible to connect the states forward, i.e. store in $S_1$ what the *next*
state is going to be: it is ambiguous!

However, when track finding has concluded, and the trajectories identified by
{term}`tip states<tip state>` $S_7$ and $S_8$ have been discarded or are being
copied into an output cotainer, it is possible to *forward link* the track
state sequences. This is possible **if** the trajectory does not branch
anymore! {func}`Acts::TrackProxy::copyFrom` will implicitly forward link the
track states, as it is guaranteed to not branch after copying.

You can manually add forward-linking to a track by calling
{func}`Acts::TrackProxy::linkForward` or
{func}`Acts::TrackProxy::reverseTrackStates`.

:::{warning}
Calling either {func}`Acts::TrackProxy::linkForward` or
{func}`Acts::TrackProxy::reverseTrackStates` on a track state sequence which
has branching will break the branching! If you have other tracks pointing at a
{term}`tip state` that branches from the sequence you're trying to
forward-link, it will be corrupted!
:::

In this example, before any forward linking, the sequence looks like this:

:::{graphviz}

digraph {
rankdir="LR";
S2 -> S1;
S3 -> S1;
S7 -> S5 -> S4 -> S2;
S6 -> S3;
}

:::

After a copy operation of $S_6$ and $S_7$ the resultig track state sequences will look like this:


:::{graphviz}

digraph {
rankdir="LR";
S11[label="S1 (copy)"];
S12[label="S1 (copy)"];
S7 -> S5 -> S4 -> S11;
S11 -> S4 -> S5 -> S7;
S6 -> S3 -> S12;
S12 -> S3 -> S6;
}

:::

This now includes both forward and backward links, which allows iteration from
$S_1$/$S_2$ to $S_6$/$S_7$ and the othe other way around.

Forward iteration can then be achieved like this:

```cpp
auto track = getTrackFromSomewhere();
for(const auto trackState : track.trackStates()) {
  // iterate forward
  // do something with track state
}
```

and the innermost track state becomes directly accessible via
{func}`Acts::TrackProxy::innermostTrackState`.

## Track State API

:::{doxygenclass} Acts::TrackStateProxy
:members:
:::

:::{doxygenclass} Acts::MultiTrajectory
:members:
:::

## Component sharing

{class}`Acts::MultiTrajectory` is designed so that components can be shared
between components. This can be achieved using the
{func}`Acts::TrackStateProxy::shareFrom` can be used to set this up.

Shareable components are
- predicted parameters and covariance
- filtered parameters and covariance
- smoothed parameters and covariance
- jacobian

To illustrate why this can be useful, consider again {numref}`ckf_tree`, where
$S_2$ and $S_3$ branch out from a shared $S_1$. In this case, the predicted
parameter vector and covariance, as well as the jacobian from $S_1\to S_2$ and
$S_1 \to S_3$ will be identical. In this case, the combinatorial track finding
will use the sharing functionality to share these components.

:::{attention}
Sharing these components introduces *cross-talk* between track states, and this
is intentional.  If e.g. the predicted covariance is modified through either of
the track states, the changes will be visible when accessed from the other
track state as well.
:::

## Dynamic columns

Aside from the static properties that both the track states and the track, the
EDM supports adding almost arbitrary additional information as dynamic columns.
The implementation of the dynamic column mechanism is given by the backend,
where the interface layer classes {class}`Acts::MultiTrajectory` and
{class}`Acts::TrackContainer` and associated proxies only coordinate the
creation, access and copying of dynamic columns.

The following illustrates the
usage of dynamic columns for {class}`Acts::TrackContainer`, but usage on
{class}`Acts::MultiTrajectory` is identical.

Assume you create a track container using some combination of backends (see
[](#edm_track_backends) for information on the backends shipped with ACTS).

```cpp
Acts::TrackContainer tc{/*...*/};

// add dynamic columns programmatically

tc.addColumn<float>("col_a");
tc.addColumn<uint8_t>("col_b");

```

Adding columns is only supported on *mutable* track containers, const track
containers should contain the original dynamic columns from when they were
created. It is up to the backend to implement recovering dynamic columns from
e.g. input files.

:::{note}
Which types are supported depends on the backend being used. See
[](#edm_track_backends) for information on the backends shipped with ACTS, and
which types they support.
:::

With these dynamic columns registered, it is now possible to set and get values for these columns on tracks.

```cpp
using namespace Acts::HashedStringLiterals;
auto track = tc.getTrack(tc.addTrack());

// these two are equivalent
track.component<float, "col_a"_hash>() = 42.42;
track.component<float>("col_a"_hash) = 52.52;
std::cout << track.component<float, "col_a"_hash>() << std::endl; // prints: 52.52
```

:::{tip}
The expression `"col_a"_hash` is a user-defined literal that internally calls
```cpp
Acts::hashedString("col_a");
```
This literal is only available after
```cpp
using namespace Acts::HashedStringLiterals;
```
:::

The components are accessed by a hash of the name of the component. This hash
can be calculated from a string at compile-time, if the string is known at
compile time.  The difference between the two component acccess signatures is
that in the first case, the hash of the component is guaranteed to be evaluated
at compile-time, since it is given to the `component` function as a template
argument. A third option is available to access components: see
[](#edm_track_accessors).


(edm_track_accessors)=
## Accessors

It can be inconvenient to have to write the full component access signature,
especially if you want to access the same components repeatedly. An alternative
are **accessors**.  They encapsulate the type of the component, and the
component name hash into an object:

```cpp
// definition of the accessor with a type and the name of the component
Acts::ProxyAccessor<float> extra("extra");
// component access by calling it on a proxy
extra(track) = 42.2;
std::cout << extra(track) << std::endl; // prints 42.2
```

:::{tip}
The same accessor also works for {class}`Acts::TrackStateProxy` objects, as it shares
the same component access mechanism with {class}`Acts::TrackProxy`.
:::

The above accessor is a **mutable** accessor, meaning it can only be used with
mutable proxy objects!

```cpp
ConstTrackProxy<...> constTrack = /*...*/;
extra(constTrack); // this will give a compile error!
```

To access properties on const proxy objects, you need to use a dedicated
accessor type:

```cpp
Acts::ConstProxyAccessor<float> extraConst("extra");
std::cout << extraConst(constTrack) << std::endl; // prints 42.2
// using the const accessor on a mutable proxy also works
std::cout << extraConst(track) << std::endl; // prints 42.2
```

For both const and mutable proxy accessors you do not actually need a mutable reference, as the internal accessor state is not mutated after construction. You can safely use a static instance of these accessors to avoid constructing them over and over again:

```cpp
template<typename track_proxy_t>
void doSomething(track_proxy_t track, float value) {
    // only created once, never changed
    static const Acts::ProxyAccessor<float> extra("extra");
    extra(track) = value;
}
```


### How to create a track from scratch

## Track EDM backends

(edm_track_backends)=
### Backends shipped with ACTS

#### Transient vector backend

:::{note}
There are currently no restrictions on types that can be used as dynamic
columns. Any type can be stored and retrieved back from the backend.

Keep in mind that the transient vector backend does not support persistency,
meaning that there is no mechanism to serialize dynamic columns (or static
columns for that matter) to disk and read them back in.
:::

#### PODIO backend

:::{note}
The PODIO backend currently supports all types that can be written to a
`podio::UserDataCollection` for dynamic columns. At the time of writing these
are: `float`, `double`, `int8_t`, `int16_t`, `int32_t`, `int64_t`, `uint8_t`,
`uint16_t`, `uint32_t`, `uint64_t`.

In particular, it is not possible to write `bool` values directly. A workaround
is using `uint8_t` and making boolean expressions explicit.

There is ongoing work to add support for PODIO *components* as values of
dynamic columns.
:::

### How to build a backend

## Glossary

:::{glossary}
tip state
  The state at the *tip* of a trajectory, usually the *outermost* track state, on the opposite end of the {term}`stem state`.

stem state
  The state at the *stem* of a trajectory, meaning the *innermost* track state. The opposite end is the {term}`tip state`.
:::
