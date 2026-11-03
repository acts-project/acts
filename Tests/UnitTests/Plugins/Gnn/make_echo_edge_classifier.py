#!/usr/bin/env python3
"""Write Tests/Data/echo_edge_classifier.onnx, used by OnnxEdgeClassifierTests.

The model ignores the graph and hands back the node features it was given,
flattened row by row, as its (numEdges x 1) edge scores. Fed numNodes *
numFeatures edges, its scores are then exactly the node feature matrix the
classifier stage passed it, which shows which features reached the model, in
which order and with which scale.

The edge index is still used, for the shape of the output, so that onnxruntime
has no reason to drop the input.
"""

from pathlib import Path

import onnx
from onnx import TensorProto, helper

nodes = [
    helper.make_node("Shape", ["edge_index"], ["edge_index_shape"]),
    helper.make_node("Slice", ["edge_index_shape", "one", "two"], ["num_edges"]),
    helper.make_node("Concat", ["num_edges", "one"], ["score_shape"], axis=0),
    helper.make_node("Reshape", ["x", "score_shape"], ["scores"]),
]
initializers = [
    helper.make_tensor("one", TensorProto.INT64, [1], [1]),
    helper.make_tensor("two", TensorProto.INT64, [1], [2]),
]
graph = helper.make_graph(
    nodes,
    "echo_edge_classifier",
    [
        helper.make_tensor_value_info("x", TensorProto.FLOAT, ["nodes", "features"]),
        helper.make_tensor_value_info("edge_index", TensorProto.INT64, [2, "edges"]),
    ],
    [helper.make_tensor_value_info("scores", TensorProto.FLOAT, ["edges", 1])],
    initializers,
)
# An old opset and IR version, so that any supported onnxruntime reads it
model = helper.make_model(graph, opset_imports=[helper.make_opsetid("", 13)])
model.ir_version = 7
onnx.checker.check_model(model)
onnx.save(
    model,
    Path(__file__).resolve().parents[3] / "Data" / "echo_edge_classifier.onnx",
)
