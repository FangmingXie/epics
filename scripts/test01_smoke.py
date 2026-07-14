"""Tier 1 smoke test: verify the SCENIC+ stack is functional in the `epics` env.

Imports the full stack and exercises the compiled/JIT core (arboreto.grnboost2 over
dask/distributed) on a tiny synthetic matrix. No real data touched.

Note: grnboost2 starts a local dask cluster that *spawns* worker processes, so the
run code MUST live under `if __name__ == '__main__':` — otherwise each spawned worker
re-imports this module and recursively launches more processes.
"""

import numpy as np
import pandas as pd

# --- import full stack (fail-fast if anything is missing) ---
import scenicplus
import pycisTopic
import pycistarget
import ctxcore
import scanpy
import anndata
import mudata
import numba
import dask

from arboreto.algo import grnboost2


def main():
    print("versions:",
          "scenicplus", scenicplus.__version__,
          "| pycisTopic", pycisTopic.__version__,
          "| numba", numba.__version__)

    N_CELLS, N_GENES, N_TFS = 50, 20, 3
    RNG = np.random.default_rng(666)

    GENE_NAMES = [f"g{i}" for i in range(N_GENES)]
    TF_NAMES = GENE_NAMES[:N_TFS]
    EXPR = pd.DataFrame(
        RNG.poisson(5.0, size=(N_CELLS, N_GENES)).astype(float),
        columns=GENE_NAMES,
    )

    adj = grnboost2(expression_data=EXPR, tf_names=TF_NAMES, seed=666, verbose=False)

    assert len(adj) > 0, "grnboost2 returned an empty adjacency table"
    assert {"TF", "target", "importance"}.issubset(adj.columns), adj.columns.tolist()
    print("grnboost2 adjacency shape:", adj.shape)
    print("SMOKE OK")


if __name__ == "__main__":
    main()
