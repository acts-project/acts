import glob
import os
import math

import pandas as pd
import numpy as np

import torch.utils

from sklearn.cluster import DBSCAN, KMeans

from sklearn.preprocessing import LabelEncoder, OrdinalEncoder
from seed_solver_network import prepareDataSet, DuplicateClassifier, Normalise


def readDataSet(CKS_files: list[str]) -> pd.DataFrame:
    """Read the dataset from the different files, remove the pure duplicate tracks and combine the datasets"""
    """
    @param[in] CKS_files: DataFrame contain the data from each track files (1 file per events usually)
    @return: combined DataFrame containing all the track, ordered by events and then by truth particle ID in each events 
    """
    data = []
    for f in CKS_files:
        datafile = pd.read_csv(f)
        datafile = prepareDataSet(datafile)
        data.append(datafile)
    return data


def prepareInferenceData(data: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Prepare the data"""
    """
    @param[in] data: input DataFrame to be prepared
    @return: array of the network input and the corresponding truth  
    """
    # Remove truth and useless variable
    target_column = "good/duplicate/fake"
    # Separate the truth from the input variables

    y = LabelEncoder().fit(data[target_column]).transform(data[target_column])
    input = data.drop(
        columns=[
            target_column,
            "seed_id",
            "Hits_ID",
            "cluster",
        ]
    )
    # Prepare the input feature
    x_cat = OrdinalEncoder().fit_transform(input.select_dtypes("object"))
    x = np.concatenate((x_cat, input), axis=1)
    return x, y


def clusterSeed(
    event: pd.DataFrame, DBSCAN_eps: float = 0.1, DBSCAN_min_samples: int = 2
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
    trackDir = event[["eta", "phi", "vertexZ", "pT"]].to_numpy()
    trackDir[:, 2] = trackDir[:, 2] / 50
    # Perform the subclustering
    clustering = DBSCAN(eps=DBSCAN_eps, min_samples=DBSCAN_min_samples).fit(trackDir)
    clusterarray = renameCluster(clustering.labels_)
    event["cluster"] = clusterarray
    sorted = event.sort_values(["cluster"], ascending=True)
    return sorted


def renameCluster(clusterarray: np.ndarray) -> np.ndarray:
    """Rename the cluster IDs to be int starting from 0"""
    """
    @param[in] clusterarray: numpy array containing the hits IDs and the cluster ID
    @return: numpy array with updated cluster IDs
    """
    last_id = -1
    new_id = -1
    for i, cluster in enumerate(clusterarray):
        if cluster != last_id or cluster == -1:
            last_id = cluster
            new_id = new_id + 1
        clusterarray[i] = new_id
    return clusterarray


# ==================================================================

import time

start = time.time()

# ttbar events as test input
CKF_files = sorted(glob.glob("odd_output" + "/event0000000[0-1][0-9]-seed.csv"))
data = readDataSet(CKF_files)

# Data of each events after clustering
clusteredData = []
# data of each events after ambiguity resolution
cleanedData = []

t1 = time.time()

# Cluster togather tracks belonging to the same particle
for event in data:
    clustered = clusterSeed(event)
    clusteredData.append(clustered)

t2 = time.time()

duplicateClassifier = torch.load("seedduplicateClassifier.pt")

import matplotlib.pyplot as plt

# Make a copy of the data to be plotted
plotData = []
for event in clusteredData:
    plotData.append(event.copy())

# Plot the distribution of the 4 variable used in the clustering
for event in plotData:
    event["eta"].hist(bins=100)
    plt.xlabel("eta")
    plt.ylabel("nb seed")
    plt.savefig("eta.png")
    plt.clf()

    event["phi"].hist(bins=100)
    plt.xlabel("phi")
    plt.ylabel("nb seed")
    plt.savefig("phi.png")
    plt.clf()

    event["vertexZ"].hist(bins=100)
    plt.xlabel("vertexZ")
    plt.ylabel("nb seed")
    plt.savefig("vertexZ.png")
    plt.clf()

    event["pT"].hist(bins=100, range=[0, 10])
    plt.xlabel("pT")
    plt.ylabel("nb seed")
    plt.savefig("pT.png")
    plt.clf()


# Create historgram filled with the number of seed per cluster
for event in plotData:
    event["nb_seed"] = 0
    event["nb_fake"] = 0
    event["nb_duplicate"] = 0
    event["nb_good"] = 0
    event["nb_cluster"] = 0
    event["nb_truth"] = 0
    event["nb_big_good"] = 0
    event["particleId"] = event.index

    event["nb_seed"] = event.groupby(["cluster"])["cluster"].transform("size")
    event["nb_seed"].hist(bins=20, weights=1 / event["nb_seed"], range=[0, 20])
    plt.xlabel("nb seed")
    plt.ylabel("nb cluster")
    plt.savefig("nb_seed.png")
    plt.clf()
    # Create historgram filled with the number of fake seed per cluster
    event.loc[event["good/duplicate/fake"] == "fake", "nb_fake"] = (
        event.loc[event["good/duplicate/fake"] == "fake"]
        .groupby(["cluster"])["cluster"]
        .transform("size")
    )
    event["nb_fake"].hist(bins=10, weights=1 / event["nb_seed"], range=[0, 10])
    plt.xlabel("nb fake")
    plt.ylabel("nb cluster")
    plt.savefig("nb_fake.png")
    plt.clf()
    # Create historgram filled with the number of duplicate seed per cluster
    event.loc[event["good/duplicate/fake"] == "duplicate", "nb_duplicate"] = (
        event.loc[event["good/duplicate/fake"] == "duplicate"]
        .groupby(["cluster"])["cluster"]
        .transform("size")
    )
    event["nb_duplicate"].hist(bins=10, weights=1 / event["nb_seed"], range=[0, 10])
    plt.xlabel("nb duplicate")
    plt.ylabel("nb cluster")
    plt.savefig("nb_duplicate.png")
    plt.clf()
    # Create historgram filled with the number of good seed per cluster
    event.loc[event["good/duplicate/fake"] == "good", "nb_good"] = (
        event.loc[event["good/duplicate/fake"] == "good"]
        .groupby(["cluster"])["cluster"]
        .transform("size")
    )
    event["nb_good"].hist(bins=10, weights=1 / event["nb_seed"], range=[0, 10])
    plt.xlabel("nb good")
    plt.ylabel("nb cluster")
    plt.savefig("nb_good.png")
    plt.clf()
    # Create historgram filled with the number of truth particle per cluster
    event["nb_truth"] = event.groupby(["cluster"])["particleId"].transform("nunique")
    event["nb_truth"].hist(bins=10, range=[0, 10])
    plt.xlabel("nb truth")
    plt.ylabel("nb cluster")
    plt.savefig("nb_truth.png")
    plt.clf()
    # Create historgram filled with the number of cluser per truth particle
    event["nb_cluster"] = event.groupby(event.index)["cluster"].transform("nunique")
    event["nb_cluster"].hist(bins=30, weights=1 / event["nb_seed"], range=[0, 30])
    plt.xlabel("nb cluster")
    plt.ylabel("nb truth")
    plt.savefig("nb_cluster.png")
    plt.clf()

    # Create historgram filled with the number of good cluser with more than one
    event["nb_good"].hist(bins=10, weights=(event["nb_seed"] > 1) / event["nb_seed"])
    plt.xlabel("nb good cluster with more than 1 seed")
    plt.ylabel("nb cluster")
    plt.savefig("nb_big_good.png")
    plt.clf()


t3 = time.time()

# Performed the MLP based ambiguity resolution
for clusteredEvent in clusteredData:
    # Prepare the data
    x_test, y_test = prepareInferenceData(clusteredEvent)
    x = torch.tensor(x_test, dtype=torch.float32)
    output_predict = duplicateClassifier(x).detach().numpy()

    # creat an array of random value between 0 and 1 of the same size as the output
    # output_predict = np.random.rand(len(x_test))

    clusteredEvent["score"] = output_predict
    # Keep only the track in cluster of more than 1 track or with a score above 0.5
    idx = (clusteredEvent["score"] > 0.0) | (
        clusteredEvent.groupby(["cluster"])["cluster"].transform("size") > 3
    )
    cleanedEvent = clusteredEvent[idx]

    # For each cluster only keep the track with the higest score
    idx = (
        cleanedEvent.groupby(["cluster"])["score"].transform(max)
        == cleanedEvent["score"]
    )
    cleanedEvent = cleanedEvent[idx]
    # cleanedEvent = cleanedEvent
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

tend = time.time()

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
print("Load : ", (t1 - start) * 1000 / len(CKF_files), "ms")
print("Clustering : ", (t2 - t1) * 1000 / len(CKF_files), "ms")
print("Inference : ", (t4 - t3) * 1000 / len(CKF_files), "ms")
print("Perf : ", (tend - t4) * 1000 / len(CKF_files), "ms")
print("tot : ", (t4 - start) * 1000 / len(CKF_files), "ms")
print("Seed filter : ", (t4 - t1) * 1000 / len(CKF_files), "ms")


# ==================================================================
# Plotting

# Combine the events together to have a better statistics
clusteredDataPlots = pd.concat(clusteredData)

cleanedDataPlots = pd.concat(cleanedData)
# cleanedDataPlots = cleanedData[0]

import matplotlib.pyplot as plt

# Plot the average score distribution for each type of track

plt.figure()
plt.hist(
    [
        cleanedDataPlots.loc[cleanedDataPlots["good/duplicate/fake"] == "good"][
            "score"
        ],
        cleanedDataPlots.loc[cleanedDataPlots["good/duplicate/fake"] == "duplicate"][
            "score"
        ],
        cleanedDataPlots.loc[cleanedDataPlots["good/duplicate/fake"] == "fake"][
            "score"
        ],
    ],
    bins=100,
    density=True,
    stacked=False,
    label=["good", "duplicate", "fake"],
)
plt.legend()
plt.xlabel("score")
plt.ylabel("number of tracks")
plt.title("Score distribution for each type of track")
plt.savefig("score_distribution.png")

# Average value of the score for 50 eta bins
averageCleanedDataPlots = cleanedDataPlots.loc[
    cleanedDataPlots["good/duplicate/fake"] == "good"
].groupby(
    pd.cut(
        cleanedDataPlots.loc[cleanedDataPlots["good/duplicate/fake"] == "good"]["eta"],
        np.linspace(-3, 3, 100),
    )
)
plt.figure()
plt.plot(
    np.linspace(-3, 3, 99),
    averageCleanedDataPlots["score"].mean(),
    label="average score",
)
plt.legend()
plt.xlabel("eta")
plt.ylabel("score")
plt.title("Average score for each eta bin")
plt.savefig("score_eta.png")

# Plot the pT distribution for each type of track
plt.figure()
plt.hist(
    [
        clusteredDataPlots.loc[clusteredDataPlots["good/duplicate/fake"] == "good"][
            "pT"
        ],
        clusteredDataPlots.loc[
            clusteredDataPlots["good/duplicate/fake"] == "duplicate"
        ]["pT"],
        clusteredDataPlots.loc[clusteredDataPlots["good/duplicate/fake"] == "fake"][
            "pT"
        ],
    ],
    bins=100,
    range=(0, 100),
    stacked=False,
    label=["good", "duplicate", "fake"],
)
plt.legend()
plt.xlabel("pT")
plt.ylabel("number of tracks")
plt.yscale("log")
plt.title("pT distribution for each type of track")
plt.savefig("pT_distribution.png")

# Average value of the score for 50 pt bins
averageCleanedDataPlots = cleanedDataPlots.loc[
    cleanedDataPlots["good/duplicate/fake"] == "good"
].groupby(
    pd.cut(
        cleanedDataPlots.loc[cleanedDataPlots["good/duplicate/fake"] == "good"]["pT"],
        np.linspace(0, 100, 50),
    )
)
plt.figure()
plt.plot(
    np.linspace(0, 100, 49),
    averageCleanedDataPlots["score"].mean(),
    label="average score",
)
plt.legend()
plt.xlabel("pT")
plt.ylabel("score")
plt.title("Average score for each eta bin")
plt.savefig("score_pt.png")
