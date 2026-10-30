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
