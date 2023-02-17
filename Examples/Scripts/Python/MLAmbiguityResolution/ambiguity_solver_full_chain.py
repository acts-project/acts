import glob
import os
import math

import pandas as pd
import numpy as np

from sklearn.preprocessing import LabelEncoder, StandardScaler, OrdinalEncoder
from sklearn.cluster import DBSCAN

import torch.nn as nn
import torch.nn.functional as F
import torch.utils

import ast
import onnxruntime as ort


def PrepareDataSet(CKS_files: list[str]) -> pd.DataFrame:
    """Read the dataset from the different file, remove the pure duplicate tracks and combine the datasets"""
    """
    @param[in] CKS_files: DataFrame contain the data from each track files (1 file per events usually)
    @return: combined DataFrame containing all the track, ordered by events and then by truth particle ID in each event 
    """
    globalindex = 0
    data = []
    for f in CKS_files:
        datafile = pd.read_csv(f)
        # Remove tracks with less than 7 measurements
        datafile = datafile[datafile["nMeasurements"] > 6]
        datafile = datafile.sort_values("good/duplicate/fake", ascending=False)
        # Remove pure duplicate (tracks purely identical) keep the ones good one if among them.
        datafile = datafile.drop_duplicates(
            subset=[
                "particleId",
                "Hits_ID",
                "nOutliers",
                "nHoles",
                "nSharedHits",
                "chi2",
            ],
            keep="first",
        )
        datafile = datafile.sort_values("particleId")
        # Set truth particle ID as index
        datafile = datafile.set_index("particleId")
        # Transform the hit list from a string to an actual list
        hitsIds = []
        mergedIds = []
        for list in datafile["Hits_ID"].values:
            hitsIds.append(ast.literal_eval(list))
        datafile["Hits_ID"] = hitsIds
        # Combine dataset
        data.append(datafile)
    return data


def PrepareData(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Prepare the data"""
    """
    @param[in] data: input DataFrame to be prepared
    @return: array of the network input and the corresponding truth  
    """
    target_column = "good/duplicate/fake"
    global avg_mean
    global avg_sdv
    # Separate the truth from the input variables
    y = LabelEncoder().fit(data[target_column]).transform(data[target_column])
    # Remove truth and useless variable
    input = data.drop(
        columns=[
            target_column,
            "track_id",
            "nMajorityHits",
            "nSharedHits",
            "truthMatchProbability",
            "Hits_ID",
            "chi2",
            "pT",
            "cluster",
        ]
    )
    # Normalise variables
    x_cat = OrdinalEncoder().fit_transform(input.select_dtypes("object"))
    x = np.concatenate((x_cat, input), axis=1)
    return x, y


def ClusterTracks(
    event: pd.DataFrame, DBSCAN_eps: float = 0.07, DBSCAN_min_samples: int = 2
) -> pd.DataFrame:
    """
    Cluster together all the track that appear to belong to the same truth particle
    To cluster the tracks together, a DBSCAN is first used followed by a sub clustering based on hits shared by tracks.
    """
    """
    @param[in] event: input DataFrame that contain all track in one event
    @param[in] DBSCAN_eps: minimum radius used by the DBSCAN to cluster track together
    @param[in] DBSCAN_min_samples: minimum number of tracks needed for DBSCAN to create a cluster
    @return: DataFrame identical to the output with an added collumn with the cluster 
    """
    # Perform the DBSCAN clustering and sort the Db by cluster ID
    trackDir = event[["eta", "phi"]].to_numpy()
    clustering = DBSCAN(eps=DBSCAN_eps, min_samples=DBSCAN_min_samples).fit(trackDir)
    # Set "cluster" to 0 if you want to see the performance without DBSCAN
    event["cluster"] = clustering.labels_
    # event["cluster"] = 0
    sorted = event.sort_values(["cluster", "nMeasurements"], ascending=[True, False])
    updatedCluster = []
    cluster_hits = sorted.loc[:, ("Hits_ID", "cluster")]
    # Further split each cluster into subCluster that have shared hits' IDs
    for key, frame in cluster_hits.groupby("cluster"):
        clusterarray = frame.to_numpy()
        clusterarray = subClustering(clusterarray, key, key)
        updatedCluster.extend(clusterarray[:, 1])
    sorted.loc[:, ("cluster")] = updatedCluster
    # Turn back the cluster ID into int
    sorted = sorted.sort_values("cluster")
    clusterarray = sorted.loc[:, ("Hits_ID", "cluster")].to_numpy()
    clusterarray = renameCluster(clusterarray)
    sorted.loc[:, ("cluster")] = clusterarray[:, 1]
    return sorted


def subClustering(clusterarray: np.ndarray, c: int, lastCluster: float) -> np.ndarray:
    """SubClustering algorithm, cluster together tracks that share hits (TODO : doesn't handle real shared hits)"""
    """
    @param[in] clusterarray: numpy array containing the hits IDs and the cluster ID
    @param[in] c: ID of the cluster we are working on
    @param[in] lastCluster: ID given to the last subcluster
    @return: numpy array with updated cluster IDs
    """
    # New cluster ID set to the float increment
    newCluster = math.nextafter(lastCluster, c + 1)
    if newCluster >= c + 1:
        raise RuntimeError(
            "Too many subcluster in the clusters, this shouldn't be possible."
        )
    hits_IDs = []
    set_IDs = set(hits_IDs)
    # Update the cluster ID of all tracks sharing a hit with the first hits that haven't been updated yet
    for track in clusterarray:
        if track[1] == c:
            if hits_IDs == []:
                hits_IDs = track[0]
                set_IDs = set(hits_IDs)
            if set_IDs & set(track[0]):
                track[1] = newCluster
    # If all ID have been updated our work is done return the updated array
    if hits_IDs == []:
        return clusterarray
    else:
        # Performe a new subclusterning for the remaning tracks
        clusterarray = subClustering(clusterarray, c, newCluster)
        return clusterarray


def renameCluster(clusterarray: np.ndarray) -> np.ndarray:
    """Rename the cluster IDs to be int starting from 0"""
    """
    @param[in] clusterarray: numpy array containing the hits IDs and the cluster ID
    @return: numpy array with updated cluster IDs
    """
    last_id = -1
    new_id = -1
    for track in clusterarray:
        if track[1] != last_id:
            last_id = track[1]
            new_id = new_id + 1
        track[1] = new_id
    return clusterarray


class DuplicateClassifier(nn.Module):
    """MLP model used to separate good tracks from duplicate tracks. Return one score per track the higher one correspond to the good track."""

    def __init__(self, input_dim, n_layers):
        """Three layer MLP, 20% dropout, sigmoid activation for the last layer."""
        super(DuplicateClassifier, self).__init__()
        self.linear1 = nn.Linear(input_dim, n_layers[0])
        self.linear2 = nn.Linear(n_layers[0], n_layers[1])
        self.linear3 = nn.Linear(n_layers[1], n_layers[2])
        self.output = nn.Linear(n_layers[2], 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, z):
        z = F.relu(self.linear1(z))
        z = F.relu(self.linear2(z))
        z = F.relu(self.linear3(z))
        return self.sigmoid(self.output(z))


class Normalise(nn.Module):
    """Normalisation of the input before the MLP model."""

    def __init__(self, mean, std):
        super(Normalise, self).__init__()
        self.mean = torch.tensor(mean, dtype=torch.float32)
        self.std = torch.tensor(std, dtype=torch.float32)

    def forward(self, z):
        z = z - self.mean
        z = z / self.std
        return z


# ==================================================================

import time

start = time.time()

import sys

sys.setrecursionlimit(10**6)

# ttbar events as test input
CKF_files = sorted(glob.glob("odd_output" + "/event0000000[0-9][0-9]-CKFtracks.csv"))
data = PrepareDataSet(CKF_files)

# Data of each events after clustering
clusteredData = []
# data of each events after ambiguity resolution
cleanedData = []

t1 = time.time()

# Cluster togather tracks belonging to the same particle
for event in data:
    clustered = ClusterTracks(event)
    clusteredData.append(clustered)

t2 = time.time()

duplicateClassifier = torch.load(os.path.dirname(__file__) + "/duplicateClassifier.pt")

t3 = time.time()


# Performed the MLP based ambiguity resolution
for clusteredEvent in clusteredData:
    # Prepare the data
    x_test, y_test = PrepareData(clusteredEvent)
    # Write the network score to a list
    output_predict = []
    for x in x_test:
        x = torch.tensor(x, dtype=torch.float32)
        output_predict.append(duplicateClassifier(x).item())

    clusteredEvent["score"] = output_predict
    cleanedEvent = clusteredEvent
    # For each cluster only keep the track with the higest score
    idx = (
        cleanedEvent.groupby(["cluster"])["score"].transform(max)
        == cleanedEvent["score"]
    )
    cleanedEvent = cleanedEvent[idx]
    cleanedData.append(cleanedEvent)

t4 = time.time()

# Compute the algorithm performances
nb_part = 0
nb_track = 0
nb_fake = 0
nb_duplicate = 0

nb_good_match = 0
nb_reco_part = 0
nb_reco_fake = 0
nb_reco_duplicate = 0
nb_reco_track = 0

for clusteredEvent, cleanedEvent in zip(clusteredData, cleanedData):
    nb_part += clusteredEvent.loc[
        clusteredEvent["good/duplicate/fake"] != "fake"
    ].index.nunique()
    nb_track += clusteredEvent.shape[0]
    nb_fake += clusteredEvent.loc[
        clusteredEvent["good/duplicate/fake"] == "fake"
    ].shape[0]
    nb_duplicate += clusteredEvent.loc[
        clusteredEvent["good/duplicate/fake"] == "duplicate"
    ].shape[0]

    nb_good_match += cleanedEvent.loc[
        cleanedEvent["good/duplicate/fake"] == "good"
    ].shape[0]
    nb_reco_fake += cleanedEvent.loc[
        cleanedEvent["good/duplicate/fake"] == "fake"
    ].shape[0]
    nb_reco_duplicate += cleanedEvent.loc[
        cleanedEvent["good/duplicate/fake"] == "duplicate"
    ].shape[0]
    nb_reco_part += cleanedEvent.loc[
        cleanedEvent["good/duplicate/fake"] != "fake"
    ].index.nunique()
    nb_reco_track += cleanedEvent.shape[0]
end = time.time()

print("===Initial efficiencies===")
print("nb particles : ", nb_part)
print("nb track : ", nb_track)
print("duplicate rate: ", 100 * nb_duplicate / nb_track, " %")
print("Fake rate: ", 100 * nb_fake / nb_track, " %")

print("===computed efficiencies===")
print("nb particles : ", nb_part)
print("nb good match : ", nb_good_match)
print("nb particle reco : ", nb_reco_part)
print("nb track reco : ", nb_reco_track)
print("Efficiency (good track) : ", 100 * nb_good_match / nb_part, " %")
print("Efficiency (particle reco) : ", 100 * nb_reco_part / nb_part, " %")
print(
    "duplicate rate: ",
    100 * ((nb_good_match + nb_reco_duplicate) - nb_reco_part) / nb_reco_track,
    " %",
)
print("Fake rate: ", 100 * nb_reco_fake / nb_reco_track, " %")

print("===computed speed===")
print("Clustering : ", (t2 - t1) * 1000 / len(CKF_files), "ms")
print("Inference : ", (t4 - t3) * 1000 / len(CKF_files), "ms")
print("tot : ", (end - start) * 1000 / len(CKF_files), "ms")

if os.path.isfile("odd_efficiency.csv"):
    resultfile = pd.read_csv("odd_efficiency.csv")
    resultfile = resultfile.set_index("PU")
    resultfile = pd.concat([resultfile, result], axis=0)
    resultfile.to_csv(path_or_buf="odd_efficiency.csv")
else:
    resultfile = result
    resultfile.to_csv(path_or_buf="odd_efficiency.csv")

for file, cleanedEvent in zip(CKF_files, cleanedData):
    newFile = file[:-4] + "-Cleaned.csv"
    cleanedEvent = cleanedEvent.sort_values("track_id")
    cleanedEvent.to_csv(path_or_buf=newFile)
