from pytorch_lightning.callbacks import ModelCheckpoint

class FilteringModelCheckpoint(ModelCheckpoint):
    def __init__(self):
        super().__init__(dirpath="checkpoints/filtering", monitor="val_loss", save_top_k=2, save_last=True, mode="min")
