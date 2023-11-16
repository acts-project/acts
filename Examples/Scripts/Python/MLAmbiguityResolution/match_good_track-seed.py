import glob

import pandas as pd
import numpy as np


def matchGood(Seed_files: list[str], CKF_files: list[str]) -> pd.DataFrame:
    """Read the dataset from the different files, remove the particle with only fakes and combine the datasets"""
    """
    @param[in] Seed_files: DataFrame contain the data from each seed files (1 file per events usually)
    @return: combined DataFrame containing all the seed, ordered by events and then by truth particle ID in each events 
    """
    data_seed = pd.DataFrame()
    data_track = pd.DataFrame()
    goodSeed = pd.DataFrame()
    data = pd.DataFrame()
    # Loop over the different track files and collect the list of seed ID associated to the good tracks
    for f_ckf, f_seed in zip(CKF_files, Seed_files):
        print("reading file: ", f_ckf, f_seed)
        data_track = pd.read_csv(f_ckf)
        data_track = data_track.loc[data_track["good/duplicate/fake"] == "good"]
        goodSeed = data_track["seed_id"]

        data_seed = pd.read_csv(f_seed)
        # Add a good seed column to the seed dataset
        data_seed["goodSeed"] = data_seed["seed_id"].isin(goodSeed)

        data_seed.loc[data_seed["good/duplicate/fake"] == "good", "good/duplicate/fake"] = "duplicate"
        data_seed.loc[data_seed["goodSeed"] == True, "good/duplicate/fake"] = "good"

        cleanedData = pd.DataFrame()

        for ID in data_seed["particleId"].unique():
            if (
                data_seed.loc[data_seed["particleId"] == ID, "goodSeed"] == False
            ).all():
                data_seed.loc[data_seed["particleId"] == ID, "good/duplicate/fake"] = "fake"
            else:
                cleanedData = pd.concat([data_seed.loc[data_seed["particleId"] == ID], cleanedData])

        # Save the cleaned dataset for future use (the cleaning is time consuming)
        matched = f_seed[:-4] + "_matched.csv"
        matchedData = data_seed.sort_values("seed_id")
        matchedData = matchedData.set_index("seed_id")
        matchedData = matchedData.drop(columns=["goodSeed"])
        matchedData.to_csv(matched)
        data = pd.concat([data, matchedData])

        # Save the cleaned dataset for future use (the cleaning is time consuming)
        cleaned = f_seed[:-4] + "_cleaned.csv"
        cleanedData = cleanedData.sort_values("seed_id")
        cleanedData = cleanedData.set_index("seed_id")
        cleanedData = cleanedData.drop(columns=["goodSeed"])
        cleanedData.to_csv(cleaned)
        
    return data


# ttbar events used as the training input
seed_files = sorted(glob.glob("odd_output" + "/event*-seed.csv"))
CKF_files = sorted(glob.glob("odd_output" + "/event*-tracks_ckf.csv"))
data = matchGood(seed_files, CKF_files)
