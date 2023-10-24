import pandas as pd
import numpy as np

import torch.nn as nn
import torch.nn.functional as F
import torch.utils

import ast


def prepareDataSet(data: pd.DataFrame) -> pd.DataFrame:
    """Format the dataset that have been written from the Csv file"""
    """
    @param[in] data: input DataFrame containing 1 event
    @return: Formated DataFrame 
    """
    data = data
    data = data.sort_values("good/duplicate/fake", ascending=False)
    # Sort by particle ID
    data = data.sort_values("particleId")
    # Set truth particle ID as index
    print(data.shape[0] / data["particleId"].nunique())
    print(
        data.loc[data["good/duplicate/fake"] == "duplicate"].shape[0]
        / data["particleId"].nunique()
    )
    print(
        data.loc[data["good/duplicate/fake"] == "fake"].shape[0]
        / data["particleId"].nunique()
    )
    print(data["particleId"].nunique())
    print(data.loc[data["good/duplicate/fake"] == "good"].shape[0])
    print("=====")
    data = data.set_index("particleId")
    # Transform the hit list from a string to an actual list
    hitsIds = []
    mergedIds = []
    for list in data["Hits_ID"].values:
        hitsIds.append(ast.literal_eval(list))
    data["Hits_ID"] = hitsIds
    # Combine dataset
    return data


class DuplicateClassifier(nn.Module):
    """MLP model used to separate goods seed from duplicate seeds. Return one score per seed the higher one correspond to the good seed."""

    def __init__(self, input_dim, n_layers):
        """Four layer MLP, sigmoid activation for the last layer."""
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


class EtaAnnotation(nn.Module):
    """Perform the eta annotation of the input before the MLP model. Based on : https://cds.cern.ch/record/2856774/files/ATL-COM-PHYS-2023-283.pdf"""

    # Doesn't appear to improve the performance of the model
    # Kept as it could be useful for further studies
    def __init__(self):
        super(EtaAnnotation, self).__init__()

        self.linear = nn.Linear(11, 30)
        self.mean = [0, 0.5, 1, 1.5, 1.8, 2, 2.2, 2.5, 3, 3.5, 4]
        self.std = [0.8, 0.7, 0.6, 0.5, 0.4, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8]

    def forward(self, z):
        eta = torch.abs(z[:, 1])
        anotation = torch.zeros(eta.size(0), 11)
        for i in range(11):
            anotation[:, i] = torch.exp(
                -((eta - self.mean[i]) ** 2) / (2 * self.std[i] ** 2)
            )
        anotation = F.relu(self.linear(anotation))
        out = torch.cat((z, anotation), 1)
        return out
