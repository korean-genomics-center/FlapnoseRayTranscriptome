# %%
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import MultipleLocator

# %%
annot_report = "/BiO/Research/FlapnoserayTumorTranscriptome/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered_att_hit_cov_abov_60.tsv"
exp_mat = "/BiO/Research/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"

# %%
def read_transcript_annot_report(annot_report):
    df_annot = pd.read_csv(annot_report, sep="\t")
    
    return df_annot

def get_prot_name(df_annot_blastp):
    dict_gene_prot_name = dict()
    dict_gene_prot_id = dict()
    for key, value in zip(df_annot_blastp["GENE"], df_annot_blastp["BLASTP"]):
        list_val = value.split(",")
        prot_name = ",".join(list_val[:-3])
        prot_id = list_val[-3]
        organism = list_val[-2]
        # if organism != "DANRE":
        #     continue
        dict_gene_prot_name[key] = prot_name
        dict_gene_prot_id[key] = prot_id
        
    return dict_gene_prot_name, dict_gene_prot_id

def read_exp_mat(exp_mat):
    df_exp = pd.read_csv(exp_mat, sep="\t", index_col=0)
    
    return df_exp

# %%
df_annot = read_transcript_annot_report(annot_report)
df_annot_blastp = df_annot[~df_annot["BLASTP"].isna()]
dict_gene_prot_name, dict_gene_prot_id = get_prot_name(df_annot_blastp)

# %%
df_exp = read_exp_mat(exp_mat)
df_exp["gene_name"] = list(map(lambda x: dict_gene_prot_id.get(x, x), list(df_exp.index)))
df_exp["gene_label"] = list(map(lambda x: dict_gene_prot_name.get(x, x), list(df_exp.index))) 
df_exp["gene_label"] = list(map(lambda x: dict_gene_prot_name.get(x, x), list(df_exp.index))) 

df_exp = df_exp.dropna(subset=["gene_name"])

# %%
df_exp = df_exp[~(df_exp["gene_name"].str.startswith("TRINITY"))]
dict_prot_id_to_name = dict(zip(df_exp["gene_name"], df_exp["gene_label"]))

# %%
df_exp_for_prerank = df_exp.copy().reset_index(drop=False)
df_exp_for_prerank["log2fc"] = df_exp_for_prerank["Left_tumor_inner"].apply(lambda x: math.log2(x+1))-df_exp_for_prerank["normal"].apply(lambda x: math.log2(x+1))
# df_exp_for_prerank["log2fc"] = df_exp_for_prerank["fc"].apply(lambda x: math.log2(x+1))
df_exp_for_prerank = df_exp_for_prerank[["index", "log2fc"]]
# df_exp_for_prerank.to_csv("/BiO/Research/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.rnk", sep="\t", index=True)

# %%
df_exp = df_exp.groupby("gene_name", as_index=False).agg({
    "Left_tumor_inner": "sum",
    "normal": "sum"
})
df_exp["gene_label"] = df_exp["gene_name"].apply(lambda x : dict_prot_id_to_name.get(x, x))

#%%
# df_exp.to_csv("/BiO/Research/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.hit.eval.1e05.coverage.abov.70.gct", sep="\t", index=False)
# %%
import math

# df_exp_cnt_filt = df_exp[df_exp["normal"] > 10]
df_exp['geo_mean_expr'] = np.log2(np.sqrt(df_exp['Left_tumor_inner'] * df_exp['normal']) + 1)
df_exp["log2fc"] = np.log2(df_exp["Left_tumor_inner"] + 1) - np.log2(df_exp["normal"] + 1)
df_exp["abslog2fc"] = df_exp["log2fc"].abs()

# df_exp_cnt_filt.to_csv("/BiO/Research/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.hit.eval.1e05.coverage.abov.70.log2fc.incl.txt", sep="\t", index=False)
# %%
import re

path_annot = "/BiO/Research/FlapnoserayTumorTranscriptome/Results/7_go/diff_expressed_tumor_normal_upreg.GOseq.enriched_att_genename.tsv"
df_annot = pd.read_csv(path_annot, sep="\t")
df_annot_select = df_annot[df_annot["category"] == "GO:0009612"]
gene_ids = list(df_annot_select["gene_ids"].values)[0]
prot_ids = list(map(lambda x: x.split(":")[-1].strip(), re.split(r',(?=TRINITY_)', gene_ids)))
prot_ids

# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize=(5, 5))

# Grey: Low expression
sns.scatterplot(
    data=df_exp,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="grey",
    color="white",
    facecolor=None,
    label="Low expression"
)

# Blue: Downregulated
df_exp_neg_filt = df_exp[
    (df_exp["geo_mean_expr"] > 2) &
    (df_exp["log2fc"] < 0) 
]
sns.scatterplot(
    data=df_exp_neg_filt,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="royalblue",
    color="white",
    facecolor=None,
    label="Downregulated"
)

# Red: Upregulated
df_exp_pos_filt = df_exp[
    (df_exp["geo_mean_expr"] > 2) &
    (df_exp["log2fc"] > 0) 
]
sns.scatterplot(
    data=df_exp_pos_filt,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="firebrick",
    color="white",
    facecolor=None,
    label="Upregulated"
)

# Orange: Upregulated & GO:0009612
df_exp_pos_filt_target = df_exp_pos_filt[df_exp_pos_filt["gene_name"].isin(prot_ids)]
sns.scatterplot(
    data=df_exp_pos_filt_target,
    x='geo_mean_expr',
    y='log2fc',
    s=100,
    alpha=1,
    edgecolor="k",
    color="orange",
    zorder=3,
    label="GO:0009612"
)

ax = plt.gca()
ax.xaxis.set_major_locator(MultipleLocator(2))

plt.axhline(0, color='black', linestyle='--', lw=2)
plt.grid(axis="both", linestyle="--", lw=0.5)
plt.xlim(-0.5, 15)
plt.ylim(-2, 2)
plt.xlabel(r'log$_2$(Mean Expression)')
plt.ylabel(r'log$_2$(Fold Change)')

# Add legend inside plot
plt.legend(
    frameon=True,
    loc="upper right",
    fontsize=12,
    facecolor='white',
    edgecolor='black'
)

# Add labels to the right of each orange point
df_exp_pos_filt_target_head = df_exp_pos_filt_target.sort_values("abslog2fc", ascending=False).head(3)
for _, row in df_exp_pos_filt_target_head.iterrows():
    ax.text(
        row['geo_mean_expr'] + 0.15,  # offset to the right
        row['log2fc'] + 0.15,
        row['gene_name'],  # label with protein name (or change to 'gene_name')
        fontsize=12,
        va='center',
        ha='left',
        weight='bold'
    )
    
plt.tight_layout()
plt.savefig("/BiO/Research/FlapnoserayTumorTranscriptome/Results/7_go/MAplot.png", dpi=300)
plt.show()
plt.close()

# %%


