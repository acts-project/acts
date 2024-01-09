Versioning and public API
=========================

ACTS uses `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_ to
indicate breaking changes in its public API. A breaking change will result in a
new major version. The `Conventional Commits
<https://www.conventionalcommits.org/en/v1.0.0/>`_ convention is used for pull
requests to track the type of changes that have been merged and ensure that the
version number is increased correctly.

Since ACTS is still under active development not all visible symbols are
automatically considered part of the public API and as such fall under the
Semantic Versioning rules. The subset of symbols that are currently part of the
public API is outlined below.

Public API
----------

At the moment only the following modules of the core library in the ``Core``
directory contribute to the public API.

- ``EventData``
- ``MagneticField``
- ``Propagator``
- ``Surfaces``
- ``Vertexing``

Within these modules, only symbols defined directly in the ``Acts`` namespace
must be considered as public.

Private API
-----------

The following components are not part of the public API but are expected to
become part of it at a later stage:

- All modules of the core library that are not explicitly listed as part of
  the public API.
- All plugins in the ``Plugins`` directory.
- The Fatras library in the ``Fatras`` directory.

The following components will never become part of the public API:

- Symbols in helper namespaces e.g. in ``Acts::detail`` must be considered
  implementation details and are not part of the public API.
- Symbols in the namespace ``Acts::Experimental`` are not part of the public API.
- The examples framework and executables in the ``Examples`` directory.
- All units tests, integration tests, (micro)-benchmarks, and related code in
  the ``Tests`` directory.
