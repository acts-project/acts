Versioning and public API
=========================

ACTS uses `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_ to
indicate breaking changes in its public API. A breaking change will result in a
new major version. The `Conventional Commits
<https://www.conventionalcommits.org/en/v1.0.0/>`_ convention is used for pull
requests to track the type of changes that have been merged and ensure that the
version number is increased correctly.

Since ACTS is still under active development, we create **major** versions
that break the API relatively frequently. We try to limit breaking changes to about
once a month.

Public API
----------

The API surface of ACTS is very large. Everything at the top-level in the ``Acts``-namespace is considered public.

.. note::
    If a private symbol, e.g. in the `Acts::detail` or `Acts::Experimental` namespaces,
    is part of an otherwise public API, this should be **considered a bug**, which we
    would appreciate being reported as an issue.

There are a number of optional components in the ``Plugins/`` folder, only some of which are considered part of the public API:

- DD4hep
- Json
- Root
- FpeMonitoring

Private API
-----------

The following components are not part of the public API but are expected to
become part of it at a later stage:

- All plugins in the ``Plugins`` directory which are not listed under *Public API*.
- The Fatras library in the ``Fatras`` directory.
- The Alignment library in the ``Alignment`` directory.

The following components will never become part of the public API:

- Symbols in helper namespaces e.g. in ``Acts::detail`` must be considered
  implementation details and are not part of the public API.
- Symbols in the namespace ``Acts::Experimental`` are not part of the public API.
- The examples framework and executables in the ``Examples`` directory.
- All units tests, integration tests, (micro)-benchmarks, and related code in
  the ``Tests`` directory.
