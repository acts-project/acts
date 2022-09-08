import pandas as pd
import numpy as np

def make_transform_dict(detector):
    all_trans = detector[["cx","cy","cz"]].to_numpy()
    all_rot = detector[["rot_xu", "rot_xv", "rot_xw", "rot_yu", "rot_yv", "rot_yw", "rot_zu", "rot_zv", "rot_zw"]].to_numpy()
    all_rot = all_rot.reshape(-1,3,3)
    all_transforms = [ (t, r) for t, r in zip(all_trans, all_rot) ]

    return dict(zip(detector.geometry_id, all_transforms))




def get_hits_from_truth(truth, detector):
    hits = truth[["hit_id", "geometry_id","tx","ty","tz","tt"]]
    hits = hits.rename(columns={"tx": "x", "ty": "y", "tz": "z", "tt": "t"})

    return hits




def get_hits_from_measurements(event_file, detector):
    geoid_to_transform = make_transform_dict(detector)
    measurements = pd.read_csv(event_file + "-measurements.csv")
    
    pos_3d = np.zeros((len(measurements.index), 3))
    
    hits = pd.DataFrame()
    
    for i, meas in enumerate(measurements.itertuples()):
        trans, rot = geoid_to_transform[meas.geometry_id]

        x = np.array([-meas.local0, -meas.local1, 0.0])
        pos_3d[i] = trans + x @ rot
    
    hits = pd.DataFrame()
    hits["geometry_id"] = measurements.geometry_id
    hits["hit_id"] = np.arange(len(pos_3d), dtype=int)
    hits["x"] = pos_3d[:,0]
    hits["y"] = pos_3d[:,1]
    hits["z"] = pos_3d[:,2]
    
    return hits



def make_simhit_mask(event_file, truth):
    simhit_map = pd.read_csv(event_file + "-measurement-simhit-map.csv")
    simhit_mask = np.zeros(len(truth.index), dtype=bool)
    simhit_mask[simhit_map.hit_id] = True
    
    return simhit_mask
    


def load_hits_and_truth_as_trackml(event_file, detector, mask_simhits=True, use_measurements_as_hits=False):
    truth = pd.read_csv(event_file + "-truth.csv")
    truth["hit_id"] = np.array(truth.index)+1
    truth = truth.drop("index", 1)

    if not use_measurements_as_hits:
        hits = get_hits_from_truth(truth, detector)
        
        # Apply simhit map, because digitization does not necessarily find every hit
        if mask_simhits:
            simhit_mask = make_simhit_mask(event_file, truth)

            hits = hits[simhit_mask]
            truth = truth[simhit_mask]
    else:
        hits = get_hits_from_measurements(event_file, detector)
        
        #simhit_mask = make_simhit_mask(event_file, truth)
        #truth = truth[simhit_mask]


    # Dummy weight
    hits["weight"] = np.ones(len(hits.index)) / len(hits.index)

    # Detector information
    hits["volume_id"] = hits.geometry_id.map(detector.set_index('geometry_id')["volume_id"].to_dict())
    hits["layer_id"] = hits.geometry_id.map(detector.set_index('geometry_id')["layer_id"].to_dict())
    hits["module_id"] = hits.geometry_id.map(detector.set_index('geometry_id')["module_id"].to_dict())
        
    return hits, truth
    
