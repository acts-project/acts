from pytorch_lightning.callbacks import ModelCheckpoint

class GNNModelCheckpoint(ModelCheckpoint):
    def __init__(self):
        super().__init__(dirpath="checkpoints/gnn", monitor="val_loss", save_top_k=2, save_last=True, mode="min")
