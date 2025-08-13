# %%
import os
from pathlib import Path

import numpy as np
import pandas as pd

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
expr_mtx = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"

# %%
df_expr = pd.read_csv(expr_mtx, sep="\t").rename(columns={"Unnamed: 0": "index"})
df_expr['geo_mean_expr'] = np.log2(np.sqrt(df_expr['Left_tumor_inner'] * df_expr['normal']) + 1)
df_expr_sum_filt = df_expr[df_expr["geo_mean_expr"] > 2]
df_expr_sum_filt = df_expr_sum_filt.set_index("index")
df_expr_sum_filt.index.name = None
df_expr_sum_filt = df_expr_sum_filt.drop(columns=["geo_mean_expr"])
expr_mtx_filt = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt.lowcount.removed"
df_expr_sum_filt.to_csv(expr_mtx_filt, sep="\t", index=True)
