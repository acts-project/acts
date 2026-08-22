@defgroup gbts Graph-Based Track Seeding
@ingroup seeding
@brief Seeding by building and filtering a graph of hit doublets

> [!tip]
> This page documents @ref Acts::Experimental::GraphBasedTrackSeeder "GBTS" as
> implemented in ACTS today. For the classical triplet-based approach that most
> of ACTS uses, see @ref seeding. GBTS is an alternative seeding strategy, not a
> layer on top of it — the two are independent entry points producing the same
> @ref Acts::SeedContainer.

## Why a graph?

Classical ACTS seeding (@ref seeding) enumerates *triplets* of space points from
a binned grid and cuts on the helix they describe. The combinatorics of that
enumeration grow steeply with occupancy, and every triplet is evaluated in
isolation — a hit that is part of a long, clean chain of compatible doublets is
treated no differently from one that is not.

**Graph-Based Track Seeding** takes the opposite order. It first builds a
*graph*: nodes are space points, and a directed edge joins two space points on
connected detector layers whenever the pair passes a set of cheap two-point
cuts. Only afterwards does it look for structure in that graph — long chains of
mutually compatible edges — and turn the best chains into seeds. The expensive
per-candidate work is therefore done once per *edge*, not once per triplet, and
the chain length itself becomes a quality signal.

The workflow has four stages, each documented below:

1. @ref gbts-nodes — sort the space points into eta/phi-ordered graph nodes.
2. @ref gbts-graph — create edges between compatible node pairs and link
   compatible edges to each other.
3. @ref gbts-cca — propagate a "level" through the edge graph to find the
   longest chains.
4. @ref gbts-extraction — follow the best chains through a Kalman-like filter
   and emit seeds.

## Geometry and layer connections {#gbts-geometry}

GBTS does not use the ACTS tracking geometry directly. It works on its own
lightweight description: a flat list of `GbtsLayer` logical layers, each
subdivided into **eta bins**.

@ref Acts::Experimental::GbtsLayerDescription gives a layer its ID, its type
(barrel or endcap) and its extent. For a barrel layer `refCoord` is the radius
and the bounds are in @f$z@f$; for an endcap it is the other way round. The layer
ID encodes the subdetector: `id / 10000 == 8` marks a pixel barrel layer, and IDs
in the `12000`–`14000` range mark strip volumes.

Which layer pairs may be joined by an edge at all is *not* hard-coded. It comes
from a **connector file** parsed by
@ref Acts::Experimental::GbtsLayerConnectionMap::fromStream, listing source
(outer) and destination (inner) layer pairs. @ref Acts::Experimental::GbtsGeometry
combines the layer descriptions with that connection map and precomputes, for
every pair of connected layers, which *eta bin* pairs are geometrically
compatible with the allowed @f$z_0@f$ range. The result is a **bin group** list,
one inner bin together with all outer bins it may connect to, kept internal to
the geometry and handed to the graph builder as its iteration schedule. It is
ordered so that outer bins are processed before the inner bins that depend on
them.

> [!note]
> The connection table is trained offline rather than written by hand.
> @ref Acts::Experimental::GbtsLayerConnectionTool accumulates layer-pair
> statistics from simulated tracks; the
> `Examples/Scripts/Python/gbts_layer_connection_training.py` script drives it.

## Graph nodes {#gbts-nodes}

@ref Acts::Experimental::GbtsNodeStorage holds the graph nodes. Space points are
fed in one at a time through `insert`, which deliberately takes **plain scalars**
rather than an ACTS container:

@snippet{trimleft} include/Acts/Seeding/GbtsNodeStorage.hpp gbts insert

This is the same pattern as @ref Acts::CylindricalSpacePointGrid — an experiment
can fill the storage straight from its own space point EDM without first copying
everything into an @ref Acts::SpacePointContainer. Overloads exist for callers
that already have @f$r@f$ and @f$\phi@f$, and for an
@ref Acts::ConstSpacePointProxy together with the columns carrying the layer
index, cluster width and local @f$y@f$ position.

`insert` assigns the node to an eta bin via
`GbtsLayer::getEtaBin` and buffers it.
`finalize` then sorts each bin by @f$\phi@f$ and materialises the nodes into a
space point container ordered by eta bin and then by @f$\phi@f$, so that every
eta bin is **one contiguous range of node indices**. A node index is therefore
all that the rest of the algorithm needs to pass around.

The per-node data the graph builder reads lives in dynamic columns on that same
container. It is deliberately packed rather than split into one array per field,
because the innermost loop reads all of it together:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts node params

@f$\tau = \cot\theta@f$. The infinite defaults disable the @f$\tau@f$ cut
entirely; only the optional machine-learning lookup table narrows them (see
@ref gbts-ml). Alongside it sits the bookkeeping the graph builder writes:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts node edge info

Each eta bin carries its node range plus the @f$\phi@f$ index used by the sliding
window. The @f$\phi@f$ index duplicates entries shifted by @f$\pm 2\pi@f$ near the
wrap-around, so the window never has to handle wrapping:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts eta bin info

## Building the graph {#gbts-graph}

The builder walks the bin groups from @ref gbts-geometry. For each inner bin it
prepares one **sliding window** in @f$\phi@f$ per connected outer bin, whose
half-width grows with the radial separation of the two bins — the further apart
they are, the more a low-@f$p_T@f$ track can bend between them. It then loops
over the inner nodes, and for each one scans only the outer nodes inside the
window.

A candidate pair @f$(n_1, n_2)@f$ becomes an edge if it survives, in order:

| Cut | Meaning |
| --- | --- |
| @f$\Delta r > @f$ `minDeltaRadius` | the two hits are radially separated enough for @f$\tau@f$ to be meaningful |
| @f$\lvert\tau\rvert < @f$ `maxAbsTau` | within the detector's angular acceptance |
| @f$\tau@f$ inside both nodes' windows | per-node acceptance from the ML lookup table |
| @f$z_0@f$ inside `[minZ0, maxZ0]`, and @f$z@f$ at the outer radius inside the ROI | the pair points back to the luminous region |
| @f$\lvert\kappa\rvert@f$ below an @f$\eta@f$-dependent bound | consistent with the @f$p_T@f$ threshold |

with @f$\tau = \Delta z/\Delta r@f$, @f$z_0 = z_1 - r_1\tau@f$ and the curvature
proxy @f$\kappa = (\phi_2-\phi_1)/\Delta r@f$.

Surviving pairs are appended to a flat edge array:

@snippet{trimleft} include/Acts/Seeding/detail/GbtsGraphTypes.hpp gbts edge

The three fit parameters `p` are @f$\{\exp(-\eta),\ \kappa,\ \phi_1 + \kappa
r_1\}@f$.

Because the inner node's edges are written contiguously, the edges *incoming* to
a node form a contiguous range, recorded in that node's
`GbtsNodeEdgeInfo`. Immediately after creating an edge
@f$(n_1, n_2)@f$, the builder scans the edges incoming to @f$n_2@f$ — that is,
edges @f$(n_2, n_3)@f$ — and links the two whenever the implied triplet is
consistent: the @f$\tau@f$ ratio, the @f$\phi@f$ continuation and the curvature
difference must all agree within tolerance. For pixel-barrel triplets an optional
@ref Acts::Experimental::GraphBasedTrackSeeder "validateTriplets" step also fits a
circle through the three points and cuts on @f$d_0@f$ and @f$p_T@f$.

Each edge stores up to `kGbtsMaxEdgeNeighbours` (6) such neighbours. Every
inner node also accumulates a 16-bit @f$z_0@f$ **histogram bitmask** of its
confirmed edges; on the innermost layer that mask is reused to reject candidates
whose @f$z_0@f$ falls in an empty bin, and nodes with no connections at all are
skipped outright.

## Connected component analysis {#gbts-cca}

With the edge graph built, a cellular automaton assigns each edge a **level**:
the length of the longest chain of linked edges ending at it. All edges start at
level 1; in each iteration an edge whose level equals that of one of its
neighbours proposes an increment, and the proposals are committed at the end of
the iteration. The sweep repeats until nothing changes, or for at most 15
iterations.

The level is the chain-length signal that drives extraction: an edge at level
@f$L@f$ is the head of a chain spanning @f$L+1@f$ space points.

## Seed extraction {#gbts-extraction}

Edges whose level clears the minimum chain length become **chain heads**, sorted
by level so the longest chains are collected first. Each head is then followed
back through the graph by @ref Acts::Experimental::GbtsTrackingFilter.

The filter is a small Kalman filter over the chain. It carries a state of two
independent parts — a quadratic in the bending plane and a linear @f$z@f$ versus
@f$r@f$ model — and at each step extrapolates to the next node, forms a
@f$\chi^2@f$ residual and rejects the branch if either component exceeds its
threshold. Where an edge has several neighbours the filter *branches*, recursing
into each; the surviving branch with the best accumulated score wins. The score
rewards each accepted hit and penalises its @f$\chi^2@f$, so the filter prefers
long chains that are also clean.

The result is a set of seed candidates. These are reduced in two passes:

- **Clone removal.** Candidates are ranked by quality, and each space point is
  assigned to the best candidate claiming it. A candidate that has lost more than
  `hitShareThreshold` of its hits to better candidates is dropped.
- **Seed splitting.** Short, central candidates are checked for self-consistency
  by fitting the circle through three different hit subsets. If the three
  curvature estimates disagree by more than `maxInvRadDiff`, the candidate is
  emitted as two shorter "drop-out" seeds instead of one.

The surviving candidates are written to the output @ref Acts::SeedContainer, with
node indices translated back to the caller's own space point indices.

## Machine-learning assisted acceptance {#gbts-ml}

When `useClusterWidthCuts` is enabled, GBTS narrows the per-node @f$\tau@f$
window using a pre-trained lookup table indexed by **pixel cluster width**. The
physical basis is that the cluster a track leaves in a pixel module grows with
the track's
incidence angle, so the cluster width alone constrains @f$\cot\theta@f$ before any
pairing is attempted — which removes candidate pairs earlier and more cheaply
than a two-point cut can.

The table carries two sets of bounds per width bin: one for clusters comfortably
inside the module, and one for clusters within `moduleEdgeTolerance` of the module
edge, where the cluster may be truncated and the width therefore underestimates
the angle. Wide clusters in the pixel endcap are dropped entirely
(`maxEndcapClusterWidth`).

> [!note]
> The lookup table is only consulted for pixel barrel layers, and the ACTS
> examples framework does not currently provide cluster widths or local
> positions, so this path is exercised only by experiment-side integrations that
> supply them through `insert`.

## Configuration {#gbts-configuration}

The main knobs on @ref Acts::Experimental::GraphBasedTrackSeeder "GraphBasedTrackSeeder::Config":

| Option | Stage | Effect |
| --- | --- | --- |
| `connectorInputFile` | @ref gbts-geometry | layer connection table; defines which layer pairs may form edges |
| `etaBinWidthOverride` | @ref gbts-geometry | override the eta bin width from the connector file (default 0.2) |
| `minPt` | @ref gbts-graph | drives the curvature and @f$\phi@f$-window bounds |
| `nMaxPhiSlice` | @ref gbts-graph | sets the base @f$\phi@f$ sliding-window width |
| `minDeltaRadius`, `maxAbsTau` | @ref gbts-graph | doublet acceptance |
| `minZ0`, `maxZ0`, `doubletFilterRZ` | @ref gbts-graph | luminous-region cuts on the doublet |
| `tauRatioCut`, `cutDPhiMax`, `cutDCurvMax` | @ref gbts-graph | edge-to-edge linking tolerances |
| `useAdaptiveCuts`, `tauRatioCorr` | @ref gbts-graph | widen the @f$\tau@f$ tolerance when a layer is skipped |
| `validateTriplets`, `d0Max` | @ref gbts-graph | circle fit on pixel-barrel triplets |
| `nMaxEdges` | @ref gbts-graph | hard cap on the edge array (2M by default); exceeding it costs efficiency |
| `matchBeforeCreate` | @ref gbts-graph | require an existing compatible incoming edge before creating one |
| `hitShareThreshold` | @ref gbts-extraction | fraction of shared hits above which a candidate is a clone |
| `maxSeedSplitEta`, `maxInvRadDiff` | @ref gbts-extraction | seed splitting |
| `addTriplets`, `maxAbsEtaAddTripelts` | @ref gbts-extraction | allow shorter chains within an @f$\eta@f$ range |
| `useClusterWidthCuts`, `lutInputFile` | @ref gbts-ml | cluster-width based @f$\tau@f$ windows |
| `maxEndcapClusterWidth`, `moduleHalfLengthY`, `moduleEdgeTolerance` | @ref gbts-ml | cluster-width acceptance and module-edge handling |
| `lrtMode` | all | Large Radius Tracking: strip layers instead of pixel, looser cuts, shorter minimum chain |

@ref Acts::Experimental::GbtsTrackingFilter "GbtsTrackingFilter::Config" separately
controls the Kalman filter — measurement resolutions, per-step @f$\chi^2@f$
ceilings and the per-hit score reward.

## Implementation pointers {#gbts-implementation}

- Seeder and configuration: @ref Acts::Experimental::GraphBasedTrackSeeder.
- Node storage: @ref Acts::Experimental::GbtsNodeStorage. The graph EDM it
  holds - `GbtsNodeParams`, `GbtsNodeEdgeInfo`, `GbtsEtaBinInfo`, `GbtsEdge` -
  is internal and lives in `Acts/Seeding/detail/GbtsGraphTypes.hpp`.
- Geometry: @ref Acts::Experimental::GbtsGeometry,
  @ref Acts::Experimental::GbtsLayerConnectionMap, and the internal `GbtsLayer`.
- Chain following: @ref Acts::Experimental::GbtsTrackingFilter and its internal
  `GbtsEdgeState`.
- Region of interest: @ref Acts::Experimental::GbtsRoiDescriptor.
- Connection-table training: @ref Acts::Experimental::GbtsLayerConnectionTool.
- Examples integration: `ActsExamples::GraphBasedSeedingAlgorithm`, driven from
  `Examples/Scripts/Python/full_chain_itk_Gbts.py`.

A GPU implementation of the same algorithm, using an equivalent
struct-of-arrays layout, lives in the traccc plugin under
`Traccc/device/common/include/traccc/gbts_seeding`.
