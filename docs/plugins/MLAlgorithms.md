# Machine leaning algorithms

Acts allows you to replace some of the components of the tracking chain by machine learning solutions. For now a replacement to the ambiguity solver is available, but when others are implemented they will be explained here.  


## Onnx plugin

To be able to perform neural network models' inferences in C++ Acts uses [onnxruntime](https://onnxruntime.ai/). An interface to use it has been implemented as an Acts plugin, to use it you will need to compile Acts with the `ACTS_PLUGIN_ONNX` option. For more detail on how to export your model to onnx please see the documentation on their [website](https://onnxruntime.ai/docs/) 

### OnnxRuntimeBase

The `OnnxRuntimeBase` class implement the inference of a standard MLP via Onnx and just require a link to the `.onnx` file containing the model one wants to use. Please note that right now the implementation of the inference in Acts only work with a single input node and a single output node. The inference can be both perform in single entry mode and in batch mode. In single entry mode, the `runONNXInference` method takes a single vector as entry, each element of the vector corresponding to one of the input features of your network. In batch mode, the input is an `Eigen::Array` with the columns corresponding to the input features and rows to the different batch inputs.

## AmbiguityResolutionMLAlgorithm

The goal of the ambiguity solver is to remove duplicated and fake tracks that remains after the CKF. To perform this cleaning, this algorithm works in three steps :

- Clustering: tracks are clustered together, one cluster ~ one truth particle
- Ranking: tracks in each cluster are ranked, the best one is kept 
- Cleaning: last pass over all the remaining tracks to remove duplicate and fake (not implemented yet)

### Clustering

The clustering is implemented with the `clusterTracks` function. Its input is a multimap of pair of track ID and vector of measurement ID. The multimap uses the number of measurement associated with the tracks as a key, this is only a trick to sort the tracks by the number of measurements efficiently. Then for each track, starting with the one which has the most measurements, we check if a cluster shares a hits with the track. If not, we create a new cluster and associate all the hits of the current track with the cluster. If yes, the track is added to that cluster (note that the hits associated to the cluster doesn't change here). After looping over all the tracks, each of them should have been associated with a cluster.   

### Ranking

At this step we have multiple clusters of tracks. We use a NN to compute a score for each track, the closer the score is to 1, the more confident we are that the track is the best track (the best one associated with the truth particle). Then for each cluster, we select the track with the highest score.

### Cleaning

Finally, for the remaining tracks, there might be some fakes and duplicates. With the first tests it doesn't seem to be the case, so the current implementation stops here. But in the future if some configuration encounters too many fake/duplicate after the ranking a simple classification NN could be used to separate the good tracks from the others. 

### How to use 

The ML based Ambiguity solver comes with a pre-trained model to be used with the ODD. Using it is extremely simple, just call the `full_chain_odd.py` with the `--MLSolver` option. If you want to try this solver with another detector, a few more step are needed. First you need to train a model. As an input you will need the multitrajectory output from the CKF (run a full chain up to the CKF then use the `CSVWriter`). At least 100 ttbar events with 200 PU are needed (1000 would be ideal). You can then run `Examples/Scripts/Python/MLAmbiguityResolution/train_ambiguity_solver.py ` on our dataset, 2/3 of the data should be used for the training and 1/3 for the validation. You will receive a `duplicateClassifier.onnx` file as an output which now can be used as an input for the `AmbiguityResolutionMLAlgorithm` with your detector.

Additionally, two other Python files are present with the training. One called `ambiguity_solver_full_chain.py` this one perform the same job as the chain we presented here but in python using Csv files in input. The difference is that they also perform an additional pre-clustering using a DBScan. The other one `ambiguity_solver_perf.py` can be used to study the performance of the ambiguity solver (ML and not) and also takes Csv files as input.