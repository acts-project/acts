# Geometry (Experimental)

```{note}
From release X onwwards, the geometry model of Acts omits the notion of layers in favour of having a pure volume and surface based geometry.
 ```

## Detector volumes an portal surfaces

```{note}
All geometry classes can only be constructed as shared objected, using the `Object::makeShared` templated methods, however, during execution (e.g. for the navigation) `const Object*` is used for performance reasons.
```

A volume (class detector volume) in Acts is defined as a set of bounding surfaces that constrain a given three-dimensional extent in space. The boundaing surfaces (class portals) are ojbects that carry links between the volumes. In such a way, a navigation stream through the detector can be followed by using the intersection of the portals to guide from volume to another. As a consequence, however, this means that no free space between volumes can exist, i.e. the boundary surfaces are between volumes and can - if necessary be shared.

A volume must hold:
 * a complete set of portal surfaces

A volume can hold:
 * a collection of contained detector surfaces
 * a collection of contained detector volumes 

 The internal objects are stored with a class ``ObjectStore<T>`` which gives acces to non-const ``std::vector<shared_ptr<T>>`` during geometry creation, and access to ``std::vector<const T*>`` for access during applications.