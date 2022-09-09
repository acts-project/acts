import torch

# Local imports
from ..utils import make_mlp
from ..filter_base import LargeFilterBaseBalanced

class PyramidFilter(LargeFilterBaseBalanced):
    def __init__(self, hparams):
        super().__init__(hparams)
        """
        Initialise the Lightning Module that can scan over different filter training regimes
        """

        # Construct the MLP architecture
        self.net = make_mlp(
            (hparams["spatial_channels"] + hparams["cell_channels"] + hparams["emb_channels"]) * 2,
            [hparams["hidden"] // (2**i) for i in range(hparams["nb_layer"])] + [1],
            layer_norm=hparams["layernorm"],
            batch_norm=hparams["batchnorm"],
            output_activation=None,
            hidden_activation=hparams["hidden_activation"],
        )
            
        
    def forward(self, x, e):
        x = self.net(torch.cat([x[e[0]], x[e[1]]], dim=-1))
        return x

