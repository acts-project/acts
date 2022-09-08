#!/bin/bash

export MODULES_ROOT=$HOME/Documents/acts_project/acts/Plugins/ExaTrkX/training/LightningModules

# First copy the files to the corresponding module dirs
cp configs/processing.yaml $MODULES_ROOT/Processing/processing.yaml
cp configs/embedding.yaml $MODULES_ROOT/Embedding/embedding.yaml
cp configs/filter.yaml $MODULES_ROOT/Filter/filter.yaml
cp configs/gnn.yaml $MODULES_ROOT/GNN/gnn.yaml

export PYTHONWARNINGS="ignore"

export EXATRKX_DATA="./data"

export PROCESSING_OUTPUT=tmp/processing_output
export EMBEDDING_OUTPUT=tmp/embedding_output
export FILTER_OUTPUT=tmp/filter_output
export GNN_OUTPUT=tmp/gnn_output
export SEGMENTING_OUTPUT=tmp/segmenting_output

export PROJECT_NAME="train-ci"

traintrack "$@" configs/pipeline.yaml

# Finally remove all yaml files
find $MODULES_ROOT -name "*.yaml" -type f -print0 | xargs -0 /bin/rm -f
