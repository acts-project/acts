@defgroup detray_plugin Detray Plugin
@ingroup plugins

To convert an enclosed part of a tracking geometry, set
@ref ActsPlugins::DetrayPayloadConverter::Config::replacementWorld to a volume in
that geometry. The converter includes that volume and all its descendants,
including gap volumes. Outgoing portal segments from the replacement world to
excluded volumes become end-of-world segments in the payload; internal links
are preserved. A descendant linking outside the selected subtree is an error.
The original tracking geometry is not modified.

The configured `beampipeVolume` must belong to the selected subtree and is still
placed at detray index zero. Leave `replacementWorld` null to convert the full
geometry. Both pointers are non-owning and must remain valid during conversion.

Set @ref ActsPlugins::DetrayPayloadConverter::Config::includeWorld to `false`
to omit the selected world itself and convert only its descendants. This also
works with the input world when `replacementWorld` is null. Portal segments
targeting the omitted world become end-of-world segments; other links to excluded
volumes remain errors. The default is `true`.

Omitting the world removes its padding, surfaces and material: tracks terminate
at its children's shell rather than crossing the world's space. The selected
world must have descendants and cannot itself be the configured beampipe.
