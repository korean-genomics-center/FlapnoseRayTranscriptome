# %%
import matplotlib.pyplot as plt
import pandas as pd
import math 
import seaborn as sns
import numpy as np

# %%
annot_report = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report_att_hit_cov_abov_60_20250725.tsv"
exp_mat = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"

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
# df_exp_for_prerank.to_csv("/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.rnk", sep="\t", index=True)

# %%
df_exp = df_exp.groupby("gene_name", as_index=False).agg({
    "Left_tumor_inner": "sum",
    "normal": "sum"
})
df_exp["gene_label"] = df_exp["gene_name"].apply(lambda x : dict_prot_id_to_name.get(x, x))

#%%
# df_exp.to_csv("/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.hit.eval.1e05.coverage.abov.70.gct", sep="\t", index=False)
# %%
import math

# df_exp_cnt_filt = df_exp[df_exp["normal"] > 10]
df_exp['geo_mean_expr'] = np.log2(np.sqrt(df_exp['Left_tumor_inner'] * df_exp['normal']) + 1)
df_exp["log2fc"] = np.log2(df_exp["Left_tumor_inner"] + 1) - np.log2(df_exp["normal"] + 1)
df_exp["abslog2fc"] = df_exp["log2fc"].abs()

# df_exp_cnt_filt.to_csv("/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.hit.eval.1e05.coverage.abov.70.log2fc.incl.txt", sep="\t", index=False)
# %%
import re

path_annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Scripts/diff_expressed_tumor_normal_upreg.GOseq.enriched_att_genename.tsv"
df_annot = pd.read_csv(path_annot, sep="\t")
df_annot_select = df_annot[df_annot["category"] == "GO:0009612"]
gene_ids = list(df_annot_select["gene_ids"].values)[0]
prot_ids = list(map(lambda x: x.split(":")[-1].strip(), re.split(r',(?=TRINITY_)', gene_ids)))
prot_ids

# %%
plt.rcParams["font.size"] = 15
plt.figure(figsize=(5, 5))

sns.scatterplot(
    data=df_exp,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="k",
    color="white"
)

df_exp_neg_filt = df_exp[
    (df_exp["geo_mean_expr"] > np.log2(0)) &
    (df_exp["log2fc"] < 0) 
]

sns.scatterplot(
    data=df_exp_neg_filt,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="k",
    color="royalblue"
)

df_exp_pos_filt = df_exp[
    (df_exp["geo_mean_expr"] > np.log2(0)) &
    (df_exp["log2fc"] > 0) 
]

sns.scatterplot(
    data=df_exp_pos_filt,
    x='geo_mean_expr',
    y='log2fc',
    s=50,
    alpha=1,
    edgecolor="k",
    color="firebrick"
)

df_exp_target_filt = df_exp[df_exp["gene_name"].isin(prot_ids)]

sns.scatterplot(
    data=df_exp_target_filt,
    x='geo_mean_expr',
    y='log2fc',
    s=100,
    alpha=1,
    edgecolor="k",
    color="orange"
)


plt.axhline(-0.1, color='black', linestyle='--', lw=2)

plt.grid(axis="both", linestyle="--")
plt.xlim(-0.5, 15)
plt.ylim(-3, 3)
plt.xlabel('log2 Geometric Mean Expression')
plt.ylabel('log2 Fold Change')
plt.tight_layout()
plt.show()
plt.close()

# %%
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MultipleLocator

# plt.rcParams["font.size"] = 16
# # Select top 10 upregulated and top 10 downregulated based on log2fc
# df_exp_cnt_filt["abslog2fc"] = df_exp_cnt_filt["log2fc"].abs()
# df_plot = df_exp_cnt_filt[df_exp_cnt_filt["abslog2fc"] > 0.5]
# top_up = df_plot.sort_values("log2fc", ascending=False).head(10)
# top_down = df_plot.sort_values("log2fc", ascending=True).head(10)

# # Set colors
# top_up["Color"] = "#d73027"
# top_down["Color"] = "#4575b4"

# # Sort for plotting
# top_up = top_up.sort_values("log2fc")
# top_down = top_down.sort_values("log2fc")

# # Create subplots
# fig, axs = plt.subplots(2, 1, figsize=(10, 10), facecolor="white")

# # Add subplot labels A and B to the left side of the figure
# fig.text(0.01, 0.95, "A", fontsize=22, fontweight='bold')  # Corresponds to top subplot
# fig.text(0.01, 0.47, "B", fontsize=22, fontweight='bold')  # Corresponds to bottom subplot

# # Plot upregulated (A)
# axs[0].barh(top_up["gene_label"], top_up["log2fc"], color=top_up["Color"])
# axs[0].set_title("Up-regulated", fontsize=18)
# axs[0].set_xlabel("log2(Fold Change)")
# axs[0].set_xlim(-3, 3)
# axs[0].axvline(0, color="gray", linestyle="--")
# axs[0].grid(axis='x', linestyle='--', alpha=0.5)

# # Plot downregulated (B)
# axs[1].barh(top_down["gene_label"], top_down["log2fc"], color=top_down["Color"])
# axs[1].set_title("Down-regulated", fontsize=18)
# axs[1].set_xlabel("log2(Fold Change)")
# axs[1].set_xlim(-3, 3)
# axs[1].axvline(0, color="gray", linestyle="--")
# axs[1].grid(axis='x', linestyle='--', alpha=0.5)

# # Formatting
# for ax in axs:
#     ax.set_xlim(-3, 3)
#     ax.xaxis.set_major_locator(MultipleLocator(1))

# plt.tight_layout()
# plt.subplots_adjust(hspace=0.3)
# plt.show()
# plt.close()
