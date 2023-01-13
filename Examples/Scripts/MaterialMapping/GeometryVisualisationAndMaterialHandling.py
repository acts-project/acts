import json

import json

def dumper(obj):
    try:
        return obj.toJSON()
    except:
        return obj.__dict__

class index_info(object):

    def __init__(self, index):
        self.index = index
        self.boundaries = []
        self.layers_no_approach = []
        self.layers_with_approach = []
        self.approaches = []
        self.passiveDiscBinningR = 0
        self.passiveDiscBinningPhi = 0
        self.passiveCylinderBinningZ = 0
        self.passiveCylinderBinningPhi = 0
        self.activeBinningRorZ = 0
        self.activeBinningPhi = 0

    def __repr__(self):
        return repr((self.index, self.boundaries, self.layers_no_approach,
                     self.layers_with_approach, self.approaches,
                     self.passiveDiscBinningR, self.passiveDiscBinningPhi,
                     self.passiveCylinderBinningZ, self.passiveCylinderBinningPhi,
                     self.activeBinningRorZ, self.activeBinningPhi))

def dump_geo(filename, plot, output_folder, dump_steering, steering_file):
    f = open(filename)
    data = json.load(f)

    # First you want to read the Volume entries and save the relevant quantities
    # - names
    index_to_names = []
    for entry in data['Volumes']['entries']:
        index_to_names.append(entry['value']['NAME'])

    # In case volume names are not found check in the volume names with dummy values
    if not index_to_names:
        for entry in data['Surfaces']['entries']:
            if "volume" in entry:
                vol = entry['volume']
                if 'volume'+str(vol) not in index_to_names:
                    index_to_names.append('volume'+str(vol))


    # Stores the information to easily decide on which layers you want the material to be mapped
    steering_map = {}

    # Once you have read the volumes extensions, you read the surfaces representing layers and boundaries.
    # All these surfaces belong to a volume, they have therefore a volume index.
    # You are interested only in:
    # - surfaces with layer index (NO approach index):
    #    -- even layer indices represent "active" layers, i.e. the ones corresponding to sensitive layers
    #    -- odd event indices represent navigation layers
    # - surfaces with layer index AND approach index:
    #    -- indicate approach layers to "active" layers:
    #       e.g. it can happen that for cylinders: 1 sits in front of the active layer,
    #                                              2 sits in the middle of the layer,
    #                                              3 sits behind the layer.
    #                                              ...
    # - surfaces with boundary index (no layer index in this case):
    #    -- they indicate boundaries between volumes. You are interested in boundaries between volumes
    #       containing at least an active layer.

    index_to_extends_layers_approach_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_approach_discs = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_bounds_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_bounds_discs = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_approach_with_material = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_discs = [[] for _ in range(len(index_to_names))]

    volume_shift = [0] * len(index_to_names)

    for entry in data['Surfaces']['entries']:
        if "layer" in entry:
            extends = []
            vol = entry['volume']

            if "sensitive" in entry:
                continue

            z_shift = 0.
            if entry['value']['transform']['translation'] != None:
                z_shift = entry['value']['transform']['translation'][2]

            if entry['value']['type'] == "CylinderSurface":
                extends = [0., entry['value']['bounds']['values'][0],
                           z_shift-entry['value']['bounds']['values'][1], z_shift+entry['value']['bounds']['values'][1]]
            else:
                extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift]

            if "approach" in entry:
                extends.append(entry["approach"])
            else:
                extends.append(-1)
            extends.append(entry["layer"])

            if entry['value']['type'] == "CylinderSurface":
                if "approach" in entry:
                    index_to_extends_layers_approach_cylinders[vol-1].append(extends)
                else:
                    index_to_extends_layers_cylinders[vol-1].append(extends)
            else:
                if "approach" in entry:
                    index_to_extends_layers_approach_discs[vol-1].append(extends)
                else:
                    index_to_extends_layers_discs[vol-1].append(extends)

        if "boundary" in entry:
            extends = []
            vol = entry['volume']

            z_shift = 0.
            if entry['value']['transform']['translation'] != None:
                z_shift = entry['value']['transform']['translation'][2]

            if entry['value']['type'] == "CylinderSurface":
                extends = [0., entry['value']['bounds']['values'][0],
                           z_shift-entry['value']['bounds']['values'][1], z_shift+entry['value']['bounds']['values'][1],entry['boundary']]
                index_to_extends_layers_bounds_cylinders[vol-1].append(extends)
            else:
                extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift,entry['boundary']]
                index_to_extends_layers_bounds_discs[vol-1].append(extends)

    # Steering the information and collect it into an output file if needed

    interesting_volumes = []
    v_index = 0
    for elements in index_to_extends_layers_cylinders:
        for coords in elements:
            if v_index not in interesting_volumes:
                interesting_volumes.append(v_index)
            if index_to_names[v_index] not in steering_map:
                steering_map[index_to_names[v_index]] = index_info(v_index+1)
            steering_map[index_to_names[v_index]].layers_no_approach.append(coords[5])
        v_index = v_index+1

    v_index = 0
    for elements in index_to_extends_layers_discs:
        for coords in elements:
            if v_index not in interesting_volumes:
                interesting_volumes.append(v_index)
            if index_to_names[v_index] not in steering_map:
                steering_map[index_to_names[v_index]] = index_info(v_index+1)
            steering_map[index_to_names[v_index]].layers_no_approach.append(coords[4])
        v_index = v_index+1

    v_index = 0
    for elements in index_to_extends_layers_bounds_cylinders:
        for coords in elements:
            if v_index in interesting_volumes:
                if index_to_names[v_index] not in steering_map:
                    steering_map[index_to_names[v_index]] = index_info(v_index+1)
                steering_map[index_to_names[v_index]].boundaries.append(coords[4])
        v_index = v_index+1

    v_index=0
    for elements in index_to_extends_layers_bounds_discs:
        for coords in elements:
            if v_index in interesting_volumes:
                if index_to_names[v_index] not in steering_map:
                    steering_map[index_to_names[v_index]] = index_info(v_index+1)
                steering_map[index_to_names[v_index]].boundaries.append(coords[3])
        v_index = v_index+1

    v_index=0
    for elements in index_to_extends_layers_approach_cylinders:
        for coords in elements:
            if coords[5] not in steering_map[index_to_names[v_index]].layers_with_approach:
                steering_map[index_to_names[v_index]].layers_with_approach.append(coords[5])
            if coords[4] not in steering_map[index_to_names[v_index]].approaches:
                steering_map[index_to_names[v_index]].approaches.append(coords[4])
        v_index = v_index+1

    v_index = 0
    for elements in index_to_extends_layers_approach_discs:
        for coords in elements:
            if coords[4] not in steering_map[index_to_names[v_index]].layers_with_approach:
                steering_map[index_to_names[v_index]].layers_with_approach.append(coords[4])
            if coords[3] not in steering_map[index_to_names[v_index]].approaches:
                steering_map[index_to_names[v_index]].approaches.append(coords[3])
        v_index = v_index+1

    if dump_steering:
        output_map = { "SteeringField" : steering_map}
        with open(steering_file, 'w', encoding='utf-8') as f:
            json.dump(output_map, f, default=dumper, ensure_ascii=False, indent=4)

    # Once you have all the data in hands, you make a few nice plots to show volumes and surfaces
    if plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'figure.max_open_warning': 0})
        from matplotlib.pyplot import cm
        from itertools import cycle
        import numpy as np

        color = cm.rainbow(np.linspace(0, 1, len(index_to_extends_layers_cylinders)))

        is_in_legend = []

        plt.figure(figsize=(20,10))
        # Plot all layers (w/o approach layers)
        v_index = 0
        for elements in index_to_extends_layers_cylinders:
            for coords in elements:
                A, B = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if index_to_names[v_index] not in is_in_legend:
                    plt.plot(x_values, y_values, c=color[v_index], label='v: '+str(v_index+1)+', '+index_to_names[v_index])
                    is_in_legend.append(index_to_names[v_index])
                else:
                    plt.plot(x_values, y_values, c=color[v_index])
            v_index = v_index+1

        v_index=0
        for elements in index_to_extends_layers_discs:
            for coords in elements:
                A, B = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if index_to_names[v_index] not in is_in_legend:
                    plt.plot(x_values, y_values, c=color[v_index], label='v: '+str(v_index+1)+', '+index_to_names[v_index])
                    is_in_legend.append(index_to_names[v_index])
                else:
                    plt.plot(x_values, y_values, c=color[v_index])
            v_index = v_index+1

        # Plot boundaries
        v_index = 0
        for elements in index_to_extends_layers_bounds_cylinders:
            for coords in elements:
                A, B = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if v_index in interesting_volumes:
                    plt.plot(x_values, y_values, c=color[v_index])
            v_index = v_index+1

        v_index=0
        for elements in index_to_extends_layers_bounds_discs:
            for coords in elements:
                A, B = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if v_index in interesting_volumes:
                    plt.plot(x_values, y_values, c=color[v_index])
            v_index = v_index+1

        plt.xlabel('z [mm]')
        plt.ylabel('R [mm]')
        plt.title("Volumes and Layers (no approach layers)")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(output_folder+'/volumes_and_layers.png')

        v_index=0
        approach_colors = ["black", "blue", "red", "green", "orange", "purple", "pink"]
        for elements in index_to_extends_layers_cylinders:
            l_index=0
            if not elements:
                v_index = v_index+1
                continue
            plt.figure(figsize=(20,10))
            color_layers = cm.rainbow(np.linspace(0, 1, len(elements)))
            for coords in elements:
                C, D = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c=color_layers[l_index], label='l: '+str(coords[5]))
                l_index=l_index+1

                for a_coords in index_to_extends_layers_approach_cylinders[v_index]:
                    if a_coords[5] == coords[5]:
                        C, D = [a_coords[3], a_coords[1]], [a_coords[2], a_coords[1]]
                        ax_values = [D[0], C[0]]
                        ay_values = [D[1], C[1]]
                        plt.plot(ax_values, ay_values, linestyle=(0, (5, 10)), c=approach_colors[a_coords[4]], label='l: '+str(coords[5])+', a: '+str(a_coords[4]))
                for a_coords in index_to_extends_layers_approach_discs[v_index]:
                    if a_coords[4] == coords[5]:
                        C,D = [a_coords[2], a_coords[0]], [a_coords[2], a_coords[1]]
                        ax_values = [D[0], C[0]]
                        ay_values = [D[1], C[1]]
                        plt.plot(ax_values, ay_values, linestyle=(0, (5, 10)), c=approach_colors[a_coords[3]], label='l: '+str(coords[5])+', a: '+str(a_coords[3]))

            plt.xlabel('z [mm]')
            plt.ylabel('R [mm]')
            plt.title(index_to_names[v_index])
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(output_folder+'/layers_for_volume_'+str(v_index+1)+'.png')
            v_index = v_index+1

        v_index=0
        for elements in index_to_extends_layers_discs:
            l_index=0
            if not elements:
                v_index = v_index+1
                continue
            plt.figure(figsize=(20,10))
            color_layers = cm.rainbow(np.linspace(0, 1, len(elements)))
            for coords in elements:
                C,D = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c=color_layers[l_index], label='l: '+str(coords[4]))
                l_index=l_index+1

                for a_coords in index_to_extends_layers_approach_discs[v_index]:
                    if a_coords[4] == coords[4]:
                        C,D = [a_coords[2], a_coords[0]], [a_coords[2], a_coords[1]]
                        ax_values = [D[0], C[0]]
                        ay_values = [D[1], C[1]]
                        plt.plot(ax_values, ay_values, linestyle=(0, (5, 10)), c=approach_colors[a_coords[3]], label='l: '+str(coords[4])+', a: '+str(a_coords[3]))
                for a_coords in index_to_extends_layers_approach_cylinders[v_index]:
                    if a_coords[5] == coords[4]:
                        C, D = [a_coords[3], a_coords[1]], [a_coords[2], a_coords[1]]
                        ax_values = [D[0], C[0]]
                        ay_values = [D[1], C[1]]
                        plt.plot(ax_values, ay_values, linestyle=(0, (5, 10)), c=approach_colors[a_coords[4]], label='l: '+str(coords[4])+', a: '+str(a_coords[4]))

            plt.xlabel('z [mm]')
            plt.ylabel('R [mm]')
            plt.title(index_to_names[v_index])
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.savefig(output_folder+'/layers_for_volume_'+str(v_index+1)+'.png')
            v_index = v_index+1

        plt.figure(figsize=(20,10))

        v_index=0
        for elements in index_to_extends_layers_bounds_cylinders:
            for coords in elements:
                A, B = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if v_index in interesting_volumes:
                    plt.plot(x_values, y_values, linestyle='--', c=color[v_index], label='v: '+str(v_index+1)+', b: '+str(coords[4]))
            v_index = v_index+1

        v_index=0
        for elements in index_to_extends_layers_bounds_discs:
            for coords in elements:
                A, B = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [B[0], A[0]]
                y_values = [B[1], A[1]]
                if v_index in interesting_volumes:
                    plt.plot(x_values, y_values, linestyle='--', c=color[v_index], label='v: '+str(v_index+1)+', b: '+str(coords[3]))
            v_index = v_index+1

        plt.xlabel('z [mm]')
        plt.ylabel('R [mm]')
        plt.title("Boundary surfaces")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(output_folder+'/boundaries.png')

        plt.figure(figsize=(20,10))

        v_index = 0
        add_to_legend=[]
        for elements in index_to_extends_layers_approach_cylinders:
            for coords in elements:
                C, D = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c=color[v_index])
                if coords[4] not in add_to_legend:
                    plt.plot(x_values, y_values, c=approach_colors[coords[4]],
                             linestyle="--", label='approach index = '+str(coords[4]))
                    add_to_legend.append(coords[4])
                else:
                    plt.plot(x_values, y_values, c=approach_colors[coords[4]],
                             linestyle="--")
            v_index = v_index+1

        v_index = 0
        for elements in index_to_extends_layers_approach_discs:
            for coords in elements:
                C,D = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c=color[v_index])
                if coords[3] not in add_to_legend:
                    plt.plot(x_values, y_values, c=approach_colors[coords[3]],
                             linestyle="--", label='approach index = '+str(coords[3]))
                    add_to_legend.append(coords[3])
                else:
                    plt.plot(x_values, y_values, c=approach_colors[coords[3]],
                             linestyle="--")
            v_index = v_index+1

        plt.xlabel('z [mm]')
        plt.ylabel('R [mm]')
        plt.title("Approach layers")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(output_folder+'/approach_layers.png')

def read_and_modify(filename, plot, output_folder, steering_file, output_file):
    f = open(filename)
    layers = open(steering_file)
    data = json.load(f)
    full_data = json.load(layers)
    layer_data = full_data["SteeringField"]

    index_to_names = []

    # Allows to dump the binning defined for the material layers
    dump_binning_for_material = False

    # Allows to check which layers are configured to carry material
    check_material_layers = True

    # First you want to read the Volume entries and save the relevant quantities
    # - names
    for entry in data['Volumes']['entries']:
        index_to_names.append(entry['value']['NAME'])

    if not index_to_names:
        index_to_names.append("FullVolume")

    index_to_extends_layers_approach_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_approach_discs = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_bounds_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_bounds_discs = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_approach_with_material = [[] for _ in range(len(index_to_names))]

    index_to_extends_layers_cylinders = [[] for _ in range(len(index_to_names))]
    index_to_extends_layers_discs = [[] for _ in range(len(index_to_names))]

    for entry in data['Surfaces']['entries']:
        if "layer" in entry:
            extends = []
            vol = entry['volume']

            z_shift = 0.
            if entry['value']['transform']['translation'] != None:
                z_shift = entry['value']['transform']['translation'][2]

            if entry['value']['type'] == "CylinderSurface":
                extends = [0., entry['value']['bounds']['values'][0],
                           z_shift-entry['value']['bounds']['values'][1], z_shift+entry['value']['bounds']['values'][1]]
            else:
                extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift]

            if "approach" in entry:
                extends.append(entry["approach"])
            else:
                extends.append(-1)
            extends.append(entry["layer"])

            if index_to_names[vol-1] in layer_data:
                if "approach" in entry:
                    if entry["layer"] in layer_data[index_to_names[vol-1]]["layers_with_approach"] and entry["approach"] in layer_data[index_to_names[vol-1]]["approaches"]:
                        entry['value']['material']['mapMaterial']=True
                else:
                    if entry["layer"] in layer_data[index_to_names[vol-1]]["layers_no_approach"]:
                        entry['value']['material']['mapMaterial']=True

                if entry['value']['material']['mapMaterial']:
                    for val in entry['value']['material']['binUtility']['binningdata']:
                        if val['value'] == 'binZ' or val['value'] == 'binR':
                            val['bins'] = layer_data[index_to_names[vol-1]]["activeBinningRorZ"]
                        else:
                            val['bins'] = layer_data[index_to_names[vol-1]]["activeBinningPhi"]
                        if val['bins']==0:
                            print("ERROR!!! Using binning value == 0! Check you input for", index_to_names[vol-1])
                            return

            if entry['value']['type'] == "CylinderSurface":
                if "approach" in entry:
                    index_to_extends_layers_approach_cylinders[vol-1].append(extends)
                    if dump_binning_for_material and entry['value']['material']['mapMaterial']:
                        print("Volume: ", entry['volume'], index_to_names[vol-1], " - Layer: ", entry['layer'], " - Approach:",entry['approach'])
                        for val in entry['value']['material']['binUtility']['binningdata']:
                            print("-->", val['value'], ": ", val['bins'])
                else:
                    index_to_extends_layers_cylinders[vol-1].append(extends)
                    if dump_binning_for_material and entry['value']['material']['mapMaterial']:
                        print("Volume: ", entry['volume'], index_to_names[vol-1], " - Layer: ", entry['layer'])
                        for val in entry['value']['material']['binUtility']['binningdata']:
                            print("-->", val['value'], ": ", val['bins'])
            else:
                if "approach" in entry:
                    index_to_extends_layers_approach_discs[vol-1].append(extends)
                    if dump_binning_for_material and entry['value']['material']['mapMaterial']:
                        print("Volume: ", entry['volume'], index_to_names[vol-1], " - Layer: ", entry['layer'], " - Approach:",entry['approach'])
                        for val in entry['value']['material']['binUtility']['binningdata']:
                            print("-->", val['value'], ": ", val['bins'])
                else:
                    index_to_extends_layers_discs[vol-1].append(extends)
                    if dump_binning_for_material and entry['value']['material']['mapMaterial']:
                        print("Volume: ", entry['volume'], index_to_names[vol-1], " - Layer: ", entry['layer'])
                        for val in entry['value']['material']['binUtility']['binningdata']:
                            print("-->", val['value'], ": ", val['bins'])

        if "boundary" in entry:
            extends = []
            vol = entry['volume']

            if index_to_names[vol-1] in layer_data and entry["boundary"] in layer_data[index_to_names[vol-1]]["boundaries"]:
                entry['value']['material']['mapMaterial']=True
                for val in entry['value']['material']['binUtility']['binningdata']:
                    if entry['value']['type'] == "CylinderSurface":
                        if val['value'] == 'binZ':
                            val['bins'] = layer_data[index_to_names[vol-1]]["passiveCylinderBinningZ"]
                        else:
                            val['bins'] = layer_data[index_to_names[vol-1]]["passiveCylinderBinningPhi"]
                    else:
                        if val['value'] == 'binR':
                            val['bins'] = layer_data[index_to_names[vol-1]]["passiveDiscBinningR"]
                        else:
                            val['bins'] = layer_data[index_to_names[vol-1]]["passiveDiscBinningPhi"]
                    if val['bins']==0:
                        print("ERROR!!! Using binning value == 0! Check you input for", index_to_names[vol-1])
                        return

            if dump_binning_for_material and entry['value']['material']['mapMaterial']:
                print("Volume: ", entry['volume'], index_to_names[vol-1], " - Boundary:",entry['boundary'])
                for val in entry['value']['material']['binUtility']['binningdata']:
                    print("-->", val['value'], ": ", val['bins'])

            z_shift = 0.
            if entry['value']['transform']['translation'] != None:
                z_shift = entry['value']['transform']['translation'][2]

            if entry['value']['type'] == "CylinderSurface":
                extends = [0., entry['value']['bounds']['values'][0],
                           z_shift-entry['value']['bounds']['values'][1], z_shift+entry['value']['bounds']['values'][1],entry['boundary']]
                index_to_extends_layers_bounds_cylinders[vol-1].append(extends)

            else:
                extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift,entry['boundary']]

                index_to_extends_layers_bounds_discs[vol-1].append(extends)

# Once you have all the data in hands, you make a few nice plots to show volumes and surfaces

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=4)

    if plot and check_material_layers:
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import cm
        from itertools import cycle
        import numpy as np

        plt.figure(figsize=(20,10))

        material_layer_cylinders = [[] for _ in range(len(index_to_names))]
        material_layer_discs = [[] for _ in range(len(index_to_names))]

        material_approach_cylinders = [[] for _ in range(len(index_to_names))]
        material_approach_discs = [[] for _ in range(len(index_to_names))]

        material_boundary_cylinders = [[] for _ in range(len(index_to_names))]
        material_boundary_discs = [[] for _ in range(len(index_to_names))]

        for entry in data['Surfaces']['entries']:

            if not entry['value']['material']['mapMaterial']:
                continue

            z_shift = 0.
            if entry['value']['transform']['translation'] != None:
                z_shift = entry['value']['transform']['translation'][2]

            if "layer" in entry:
                extends = []
                vol = entry['volume']

                if entry['value']['type'] == "CylinderSurface":
                    extends = [0., entry['value']['bounds']['values'][0],
                               z_shift-entry['value']['bounds']['values'][1],
                               z_shift+entry['value']['bounds']['values'][1]]

                    if "approach" in entry:
                        material_approach_cylinders[vol-1].append(extends)
                    else:
                        material_layer_cylinders[vol-1].append(extends)
                else:
                    extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift]

                    if "approach" in entry:
                        material_approach_discs[vol-1].append(extends)
                    else:
                        material_layer_discs[vol-1].append(extends)

            if "boundary" in entry:
                extends = []
                vol = entry['volume']

                if entry['value']['type'] == "CylinderSurface":
                    extends = [0., entry['value']['bounds']['values'][0],
                               z_shift-entry['value']['bounds']['values'][1],
                               z_shift+entry['value']['bounds']['values'][1]]

                    material_boundary_cylinders[vol-1].append(extends)

                else:
                    extends = [entry['value']['bounds']['values'][0], entry['value']['bounds']['values'][1],z_shift]
                    material_boundary_discs[vol-1].append(extends)

        v_index=0
        is_first=True
        for elements in material_layer_cylinders:
            l_index = 0
            for coords in elements:
                C, D = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                if is_first:
                    plt.plot(x_values, y_values, c="black", label="layer")
                    is_first=False
                else:
                    plt.plot(x_values, y_values, c="black")
                l_index = l_index+1
            v_index = v_index+1

        v_index=0
        for elements in material_layer_discs:
            l_index = 0
            for coords in elements:
                C,D = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c="black")
                l_index = l_index+1
            v_index = v_index+1

        v_index=0
        is_first=True
        for elements in material_approach_cylinders:
            l_index = 0
            for coords in elements:
                C, D = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                if is_first:
                    plt.plot(x_values, y_values, c="red", label="approach")
                    is_first=False
                else:
                    plt.plot(x_values, y_values, c="red")
                l_index = l_index+1
            v_index = v_index+1

        v_index=0
        for elements in material_approach_discs:
            l_index = 0
            for coords in elements:
                C,D = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c="red")
                l_index = l_index+1
            v_index = v_index+1

        v_index=0
        is_first=True
        for elements in material_boundary_cylinders:
            l_index = 0
            for coords in elements:
                C, D = [coords[3], coords[1]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                if is_first:
                    plt.plot(x_values, y_values, c="blue", label="boundary")
                    is_first=False
                else:
                    plt.plot(x_values, y_values, c="blue")
                l_index = l_index+1
            v_index = v_index+1

        v_index=0
        for elements in material_boundary_discs:
            l_index = 0
            for coords in elements:
                C,D = [coords[2], coords[0]], [coords[2], coords[1]]
                x_values = [D[0], C[0]]
                y_values = [D[1], C[1]]
                plt.plot(x_values, y_values, c="blue")
                l_index = l_index+1
            v_index = v_index+1

        plt.xlabel('z [mm]')
        plt.ylabel('R [mm]')
        plt.title("Layers with material")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(output_folder+'/material_layers.png')



import argparse
import os
# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("--geometry", default="", type=str,
                    help="Specify the input file for geometry visualisation")
parser.add_argument("--plot", default=True, action="store_true",
                    help="Enable plot creation for geometry visualisation (Default : True)")
parser.add_argument("--output_folder", default="plots", type=str,
                    help="Specify the output folder for plots (Default : plots)")
parser.add_argument("--dump_steering", default=False, action="store_true",
                    help="Enable production of steering file for material mapping (Default : False)")
parser.add_argument("--edit", default=False, action="store_true",
                    help="Enable editing of input file for creation of json for material mapping (Default : False)")
parser.add_argument("--steering_file", default="", type=str,
                    help="Specify the steering file to guide editing of json for material mapping")
parser.add_argument("--output_map", default="", type=str,
                    help="Specify the output json for material mapping")
args = parser.parse_args()

print (" --- Geometry visualisation and material map fie creation --- ")
print ()

if not args.geometry:
    print ("Error: Missing input geometry file. Please specify --geometry.")
    exit()

if not os.path.exists(args.geometry):
    print ("Error: Invalid file path/name in --geometry. Please check your input!")
    exit()

if not args.geometry.endswith(".json"):
    print ("Error: Invalid file format in --geometry. Please check your input!")
    exit()

print ("\t parsing file : ", args.geometry)
if args.plot:
    print ("\t job configured to produce plots in output folder: ", args.output_folder)
    if not os.path.exists(args.output_folder):
        os.mkdir(args.output_folder)

if args.dump_steering and args.edit:
    print ("Error: Wrong job configuration. --dump_steering and --edit can't be \
        both true at the same time.")
    print ("\t Decide if you want to dump the steering file OR to read an existing file for editing the geometry file.")

if args.dump_steering:
    if not args.steering_file:
        print ("Error: Missing output steering file. Please specify --steering_file.")
        exit()
    if not args.steering_file.endswith(".json"):
        print ("Error: Invalid file format in --steering_file. It must end with .json!")
        exit()
    print ("\t job configured to produce steering file for material mapping with name: ", args.steering_file)

if args.edit:
    print ("\t job configured to edit the input geometry file following a steering file")
    if not args.steering_file:
        print ("Error: Missing input steering file. Please specify --steering_file.")
        exit()
    if not os.path.exists(args.steering_file):
        print ("Error: Invalid file path/name in --steering_file. Please check your input!")
        exit()
    if not args.steering_file.endswith(".json"):
        print ("Error: Invalid file format in --steering_file. Please check your input!")
        exit()
    if not args.output_map:
        print ("Error: Missing output map file. Please specify --output_map.")
        exit()
    if not args.output_map.endswith(".json"):
        print ("Error: Invalid file format in --output_map. Please check your input!")
        exit()


if args.plot or args.dump_steering:
    dump_geo(args.geometry,
             args.plot, args.output_folder,
             args.dump_steering, args.steering_file)

if args.edit:
    read_and_modify(args.geometry,
                    args.plot, args.output_folder,
                    args.steering_file, args.output_map)
