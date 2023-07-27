import acts
import math
import acts.examples


from acts import (
    logging,
    Binning,
    ProtoBinning,
    LayerStructureBuilder,
    VolumeStructureBuilder,
    DetectorVolumeBuilder,
    DetectorBuilder,
    CylindricalContainerBuilder,
    Transform3,
)


class DetectorVolume:
    def __init__(
        self,
        name,
        shape,
        bounds,
        transform,
        surfaces=None,
        binnings=None,
        supports=None,
    ):
        """Creates a Detector Volume definition with the given dimensions, and optional internal structure

        Parameters
        ----------
        name: string
            The name of the volume
        shape: acts.VolumeBoundsType
            The bounds type of the volume
        bounds : list of float
            The bound parameters of the volume, according to the volume bounds type.
        transform : acts.Transform3
            The positioning of the volume in global coordinates
        surfaces : list of Surface
            The surfaces of the layer structure (if existing)
        binnings : list of int
            The binnings of the layer structure (if existing)
        supports : list of support structure parameters
            The supports of the layer structure (if existing)

        Returns
        -------
        DetectorVolume
            The created volume proxy that will provide an appropriate builder.
        """
        self.__name = name
        self.__shape = shape
        self.__bounds = bounds
        self.__transform = transform
        self.__surfaces = surfaces
        self.__binnings = binnings
        self.__supports = supports

    def print(self):
        print("LayerVolume: ", self.__name)
        print("  shape: ", self.__shape)
        print("  bounds: ", self.__bounds)
        print("  transform: ", self.__transform)
        print("  surfaces: ", self.__surfaces)
        print("  binnings: ", self.__binnings)
        print("  supports: ", self.__supports)

    def getBounds(self):
        """Returns the bounds of the layer volume.

        Returns
        -------
        list of float
            The dimensions of the layer volume.
        """
        return self.__bounds

    def setBounds(self, bounds):
        """Sets the bounds of the layer volume.

        Parameters
        ----------
        dimensions : list of float
           The dimensions of the layer volume.
        """
        self.__bounds = bounds

    def getSurfaces(self):
        """Returns the surfaces of the layer structure (if existing).

        Returns
        -------
        list of Surface
            The surfaces of the layer structure.
        """
        return self.__surfaces

    def getBinnings(self):
        """Returns the binnings of the layer structure (if existing).

        Returns
        -------
        list of int
            The binnings of the layer structure.
        """
        return self.__binnings

    def getSupports(self):
        """Returns the supports of the layer volume (if existing)..

        Returns
        -------
        list of Surface
            The supports of the layer volume.
        """
        return self.__supports

    def setTransform(self, transform):
        """Sets the transform of the layer volume.

        Parameters
        ----------
        transform : acts.Transform3
           The transform of the layer volume.
        """
        self.__transform = transform

    def getTransform(self):
        """Returns the transform of the layer volume.

        Returns
        -------
        acts.Transform3
            The transform of the layer volume.
        """
        return self.__transform

    def setName(self, name):
        self.__name = name

    def getName(self):
        return self.__name

    def getBuilder(self):
        """Returns the builder for the layer volume

        This creates the internal and external structure
        builder and combines them into a volume builder.
        """

        # Internal structure (if existing)
        layerStructureBuilder = None
        if self.__surfaces is not None:
            layerStructureConfig = acts.LayerStructureBuilder.Config()
            layerStructureConfig.surfacesProvider = self.__surfaces
            layerStructureConfig.binnings = self.__binnings

            layerStructureBuilder = acts.LayerStructureBuilder(
                layerStructureConfig,
                self.__name + "_LayerBuilder",
                acts.logging.VERBOSE,
            )

        # External structure
        volumeStructureConfig = acts.VolumeStructureBuilder.Config()
        volumeStructureConfig.boundsType = self.__shape
        volumeStructureConfig.boundValues = self.__bounds
        volumeStructureConfig.transform = self.__transform

        volumeStructureBuilder = acts.VolumeStructureBuilder(
            volumeStructureConfig, self.__name + "_ShapeBuilder", acts.logging.VERBOSE
        )

        # Volume builder
        volumeBuilderConfig = acts.DetectorVolumeBuilder.Config()
        volumeBuilderConfig.name = self.__name
        volumeBuilderConfig.internalsBuilder = layerStructureBuilder
        volumeBuilderConfig.externalsBuilder = volumeStructureBuilder

        volumeBuilder = acts.DetectorVolumeBuilder(
            volumeBuilderConfig, self.__name + "_VolumeBuilder", acts.logging.VERBOSE
        )

        return volumeBuilder


class ContainerStructure:
    def __init__(self, name, volumes, binning):
        """Creates a container from a given list of volumes

        Parameters
        ----------
        name: string
            The name of the container
        volumes: list of DetectorVolume
            The volumes of the container
        binning: acts.Binning
            The ordering of the container

        Returns
        -------
        ContainerStructure
            The created container proxy that will provide an appropriate builder.
        """
        self.__name = name
        self.__volumes = volumes
        self.__binning = binning

    def getName(self):
        return self.__name

    def getBuilder(self):
        """Returns the builder for the container"""
        volumeBuilders = []
        for volume in self.__volumes:
            volumeBuilders += [volume.getBuilder()]

        containerBuilder = None
        if self.__binning == acts.Binning.r or self.__binning == acts.Binning.z:
            containerConfig = acts.CylindricalContainerBuilder.Config()
            containerConfig.builders = volumeBuilders

            containerBuilder = acts.CylindricalContainerBuilder(
                containerConfig, self.__name + "Builder", acts.logging.VERBOSE
            )
        return containerBuilder


def createCylindricalContainer(name, dimensions, layer_volumes, binning):
    """Convenience function to create a cylindrical container structure

    Gap volumes will be created for this

    Parameters
    ----------
    name: string
        The name of the container builder
    dimensions: list of float
        The dimensions of the container, rmin, rmax, zmin, zmax
    dimensions: list of DetectorVolume
        The volumes of the container
    binning: acts.Binning
        The ordering of the container
    """

    volume_r_min, volume_r_max, volume_z_min, volume_z_max = dimensions
    complete_volumes = []

    # Consistency check and fill gaps
    last_gap = None
    last_r_max = volume_r_min
    last_z_max = volume_z_min
    igap = 0
    # Loop over layers and create gaps
    for layer in layer_volumes:
        # Gather layer volume dimensions
        layer_bounds = layer.getBounds()
        layer_r_min = layer_bounds[0]
        layer_r_max = layer_bounds[1]
        layer_half_z = layer_bounds[2]
        layer_z = layer.getTransform().getTranslation()[2]
        layer_z_min = layer_z - layer_half_z
        layer_z_max = layer_z + layer_half_z

        if layer_r_min < volume_r_min or layer_r_max > volume_r_max:
            raise Exception(
                "Layer volume r dimensions are bigger than container dimensions"
            )
        if layer_z_min < volume_z_min or layer_z_max > volume_z_max:
            raise Exception(
                "Layer volume z dimensions are bigger than container dimensions"
            )

        # The container is r-ordered
        if binning is acts.Binning.r:
            # Synchronize layer volume r dimensions with container
            volume_half_z = 0.5 * (volume_z_max - volume_z_min)
            volume_z = 0.5 * (volume_z_max + volume_z_min)
            volume_transform = acts.Transform3([0.0, 0.0, volume_z])

            # Force detector container dimension onto layer
            layer.setBounds([layer_r_min, layer_r_max, volume_half_z])
            layer.setTransform(volume_transform)

            if layer_r_min > last_r_max:
                gap_bounds = [last_r_max, layer_r_min, volume_half_z]
                gap_volume = DetectorVolume(
                    name + "_gap_" + str(igap),
                    acts.VolumeBoundsType.Cylinder,
                    gap_bounds,
                    volume_transform,
                )
                complete_volumes += [gap_volume, layer]
            else:
                complete_volumes += [layer]
            # potential last gap
            if layer_r_max < volume_r_max:
                gap_bounds = [layer_r_max, volume_r_max, volume_half_z]
                last_gap = DetectorVolume(
                    name + "_gap_" + str(igap + 1),
                    acts.VolumeBoundsType.Cylinder,
                    gap_bounds,
                    volume_transform,
                )
            else:
                last_gap = None
            # Remember the last
            last_r_max = layer_r_max
        # The container is z-ordered
        elif binning is acts.Binning.z:
            # Force detector container dimension onto layer
            layer.setBounds(
                [volume_r_min, volume_r_max, 0.5 * (layer_z_max - layer_z_min)]
            )
            if layer_z_min > last_z_max:
                gap_bounds = [
                    volume_r_min,
                    volume_r_max,
                    0.5 * (layer_z_min - last_z_max),
                ]
                gap_transform = acts.Transform3(
                    [0.0, 0.0, 0.5 * (last_z_max + layer_z_min)]
                )
                gap_volume = DetectorVolume(
                    name + "_gap_" + str(igap),
                    acts.VolumeBoundsType.Cylinder,
                    gap_bounds,
                    gap_transform,
                )
                complete_volumes += [gap_volume, layer]
            else:
                complete_volumes += [layer]
            # potential last gap
            if layer_z_max < volume_z_max:
                gap_bounds = [
                    volume_r_min,
                    volume_r_max,
                    0.5 * (volume_z_max - layer_z_max),
                ]
                gap_transform = acts.Transform3(
                    [0.0, 0.0, 0.5 * (volume_z_max + layer_z_max)]
                )
                last_gap = DetectorVolume(
                    name + "_gap_" + str(igap + 1),
                    acts.VolumeBoundsType.Cylinder,
                    gap_bounds,
                    gap_transform,
                )
            else:
                last_gap = None
            # Remember the last
            last_z_max = layer_z_max
        else:
            raise Exception("Binning not supported")
        igap += 1
    # Fill the last gap if existing and return
    if last_gap is not None:
        complete_volumes += [last_gap]
    # The container is successfully created
    return ContainerStructure(name, complete_volumes, binning)


# Convenience function to create a cylindrical barrel from a surface map
def createBarrel(name, dimensions, layerSetups, surfaceSetups, binningSetups):
    """Convenience function to create a cylindrical barrel from a surface map

    Parameters
    ----------
    name: string
        The name of the container builder
    dimensions: list of float
        The dimensions of the container, rmin, rmax, zmin, zmax
    """

    volume_r_min, volume_r_max, volume_z_min, volume_z_max = dimensions
    volume_half_z = 0.5 * (volume_z_max - volume_z_min)
    volume_transform = acts.Transform3([0.0, 0.0, 0.5 * (volume_z_min + volume_z_max)])

    # The barrel detector
    layers = []
    for i, layer_setup in enumerate(layerSetups):
        # The layer setup: r, half thickness in r, layerID
        r, rt = layer_setup
        layer_surfaces = surfaceSetups[i]
        binning_z, binning_phi = binningSetups[i]
        # Create and add the layer
        layer = DetectorVolume(
            name + "_layer_" + str(i),
            acts.VolumeBoundsType.Cylinder,
            [r - rt, r + rt, volume_half_z],
            volume_transform,
            layer_surfaces,
            [binning_z, binning_phi],
            [],
        )
        layers += [layer]

    # return a cylindrical container
    return createCylindricalContainer(name, dimensions, layers, acts.Binning.r)


def createEndcap(name, dimensions, layerSetups, surfaceSetups, binningSetups):
    """Creates a barrel structure with the given dimensions and layers."""
    endcapRmin, endcapRmax, endcapZmid, endcapZmax = dimensions

    layers = []
    for i, layer_setup in enumerate(layerSetups):
        # The layer setup : z, half_z, layerID
        z, half_z = layer_setup
        layer_surfaces = surfaceSetups[i]
        binning_r, binning_phi = binningSetups[i]
        # Create and add the layer
        layer = DetectorVolume(
            name + "_layer_" + str(i),
            acts.VolumeBoundsType.Cylinder,
            [endcapRmin, endcapRmax, half_z],
            acts.Transform3([0.0, 0.0, z]),
            layer_surfaces,
            [binning_r, binning_phi],
            [],
        )
        layers += [layer]

    # create a cylindrical container
    return createCylindricalContainer(name, dimensions, layers, acts.Binning.z)


def phiBinning(phiBins, extraBins):
    """helper method to create phi binning"""
    return acts.ProtoBinning(
        acts.Binning.phi, acts.Binning.closed, -math.pi, math.pi, phiBins, extraBins
    )
