import numpy as np
from pytorch_lightning.callbacks import Callback

def format_number_human(num):
    if num >= 1e9:
        return "{:.1f}G".format(num / 1.0e9)
    elif num >= 1e6:
        return "{:.1f}M".format(num / 1.0e6)
    elif num >= 1e3:
        return "{:.1f}K".format(num / 1.0e3)
    else:
        return "{:.1f}".format(num)



class GraphSizeMonitor(Callback):
    def __init__(self):
        super().__init__()

    def evaluate(self, data, name):
        nodes = []
        edges = []
        truth = []

        for graph in data:
            nodes.append(len(graph.x))
            edges.append(graph.edge_index.shape[1])
            truth.append(len(graph.y))

        print("Average {} graph info:".format(name))
        print("- nodes:",format_number_human(np.mean(nodes)))
        print("- edges:",format_number_human(np.mean(edges)))
        print("- truth:",format_number_human(np.mean(truth)))


    def on_fit_start(self, trainer, pl_module):
        self.evaluate(pl_module.trainset, "TRAIN")

    def on_test_start(self, trainer, pl_module):
        self.evaluate(pl_module.testset, "TEST")
