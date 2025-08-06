# %%
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# %%
def make_dict_enrich(enrich, qval_threshold):
    dict_enrich_go_bp = dict()
    dict_enrich_go_mf = dict()
    dict_enrich_go_cc = dict()
    with open(enrich, mode='r') as fr:
        _ = fr.readline()
        for line in fr:
            record = line.rstrip().split("\t")
            numIntGene = int(record[3])
            numUnivGene = int(record[4])
            if numUnivGene < 20:
                continue
            geneRatio = float(numIntGene / numUnivGene)
            goTerm = record[5]
            goCat = record[6]
            goCode = record[0]
            goName = f"{goTerm}({goCode})"
            qValue = float(record[7]) 

            if qValue < float(qval_threshold):
                if goCat == "BP":
                    dict_enrich_go_bp[goName] = [numIntGene, numUnivGene, geneRatio, qValue, -math.log10(qValue+1e-30)]
                
                if goCat == "MF":
                    dict_enrich_go_mf[goName] = [numIntGene, numUnivGene, geneRatio,  qValue, -math.log10(qValue+1e-30)]

                if goCat == "CC":
                    dict_enrich_go_cc[goName] = [numIntGene, numUnivGene, geneRatio, qValue, -math.log10(qValue+1e-30)]

    return dict_enrich_go_bp, dict_enrich_go_mf, dict_enrich_go_cc

# %%
qval_threshold = 0.0001
up_file = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Scripts/diff_expressed_tumor_normal_upreg.GOseq.enriched"
down_file = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Scripts/diff_expressed_tumor_normal_downreg.GOseq.enriched"
dict_up_bp, _, _ = make_dict_enrich(up_file, qval_threshold)
dict_down_bp, _, _ = make_dict_enrich(down_file, qval_threshold)

# %%
# Generate dataframes
def prepare_df(dict_enrich, goType, top=20):
    df = pd.DataFrame.from_dict(dict_enrich, orient="index").reset_index()
    df.columns = [goType, "numInputGene", "numUnivGene", "GeneRatio", "qValue", "-log10(qValue)"]
    # df = df[df["numUnivGene"] > 100]
    
    if top == None:
        df = df.sort_values(by=["GeneRatio", "qValue", "-log10(qValue)"], ascending=False)
    else:
        df = df.sort_values(by=["GeneRatio", "qValue", "-log10(qValue)"], ascending=False).head(top)
    
    df = df.reset_index(drop=True)
    
    return df

df_up_top = prepare_df(dict_up_bp, "BP", top=20)
df_up_all = prepare_df(dict_up_bp, "BP", top=None)
df_down_top = prepare_df(dict_down_bp, "BP", top=20)
df_down_all = prepare_df(dict_down_bp, "BP", top=None)

# Extract GO:0009612 entry
df_go_extra = df_up_all[df_up_all["BP"].str.contains("GO:0009612")].copy()
df_go_extra["BP"] = df_go_extra["BP"].values[0]

# Combine top 20 + GO:0009612
df_up_combined = pd.concat([df_up_top, df_go_extra], ignore_index=True)

# %%
def plot_lollipop(ax, df, title, color):
    df = df.sort_values("GeneRatio", ascending=True)
    
    for i, row in df.iterrows():
        ax.plot([-0.1, row["GeneRatio"]], [row["BP"], row["BP"]], 
                color=color, lw=3, zorder=1)

    sizes = df["-log10(qValue)"] * 20
    scatter = ax.scatter(x=df["GeneRatio"], y=df["BP"], 
                         s=sizes, c=color, edgecolors='black', lw=0.8, zorder=2)

    ax.set_title(title, fontsize=18)
    ax.set_xlabel("GeneRatio")
    ax.set_xlim(0, 1)
    ax.set_ylabel("")
    ax.grid(True, axis='x', linestyle='--', alpha=0.6)

    return scatter

# %%
from matplotlib.ticker import MultipleLocator

plt.rcParams["font.size"] = 13
fig, axs = plt.subplots(2, 1, figsize=(10, 12), facecolor="white")

sc_up = plot_lollipop(axs[0], df_up_combined, "Up-regulated", color="#d73027")
sc_down = plot_lollipop(axs[1], df_down_top, "Down-regulated", color="#4575b4")

d = .25
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
axs[0].plot([0, 0], [0.055, 0.075], transform=axs[0].transAxes, **kwargs, zorder=5)
axs[0].plot([0, 0], [0.055, 0.075], transform=axs[0].transAxes, color='white', linewidth=1, clip_on=False, zorder=3)

fig.text(0.02, 0.99, "A", fontsize=22, fontweight='bold', va='top', ha='left')
fig.text(0.02, 0.49, "B", fontsize=22, fontweight='bold', va='top', ha='left')

# Custom legend
legend_sizes = [1, 5, 10]
legend_dots = [
    plt.scatter([], [], s=s*20, edgecolors='black', facecolor='black', label=f"{s}", lw=0.5)
    for s in legend_sizes
]

# Add legends to each subplot
axs[0].legend(handles=legend_dots, title='-log10(qVal.)', loc='lower right', frameon=False)
axs[1].legend(handles=legend_dots, title='-log10(qVal.)', loc='lower right', frameon=False)

# Set x-axis tick spacing to 0.1
for ax in axs:
    ax.xaxis.set_major_locator(MultipleLocator(0.2))

plt.tight_layout()
plt.subplots_adjust(hspace=0.2)
# plt.savefig("")
plt.show()
plt.close()

# %%
# path_enrichment_upregulated_genes = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/6_deg/20250711/enrichment_upregulated_terms.xlsx"
# path_enrichment_downregulated_genes = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/6_deg/20250711/enrichment_downregulated_terms.xlsx"
# df_up_all.to_excel(path_enrichment_upregulated_genes, sep="\t")
# df_down_all.to_excel(path_enrichment_downregulated_genes, sep="\t")

# %%
path_obo = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/go-basic.obo"
from goatools.obo_parser import GODag
go_dag = GODag(path_obo)
go_depth = {go: go_dag[go].depth for go in go_dag}

# %%
list_GO_terms_up = df_up_all["BP"].to_list()
list_GO_terms_down = df_down_all["BP"].to_list()

list_GO_ID_up = list(map(lambda x: x.split("(")[-1][:-1], list_GO_terms_up))
list_GO_ID_down = list(map(lambda x: x.split("(")[-1][:-1], list_GO_terms_down))

dict_GO_ID_terms_up = dict(zip(list_GO_ID_up, list_GO_terms_up))
dict_GO_ID_terms_down = dict(zip(list_GO_ID_down, list_GO_terms_down))

# %%
list_go_depth_up = list()
for go_id, go_term in dict_GO_ID_terms_up.items():
    depth = go_depth.get(go_id, np.nan)
    list_go_depth_up.append(depth)

list_go_depth_up

list_go_depth_down = list()
for go_id, go_term in dict_GO_ID_terms_down.items():
    depth = go_depth.get(go_id, np.nan)
    list_go_depth_down.append(depth)

list_go_depth_down

# %%
import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 14
plt.hist(list_go_depth_up, color="red", alpha=0.5, label="upregulated", bins=12)
plt.hist(list_go_depth_down, color="blue", alpha=0.5, label="downregulated", bins=12)
plt.xticks(list(range(0, 13)))
plt.xlabel("GO hierarchy level", fontsize=15)
plt.ylabel("GO term count", fontsize=15)
plt.axvline(x=3, color="k", linestyle="--")
plt.axvline(x=6, color="k", linestyle="--")
plt.legend(bbox_to_anchor=(1, 1), frameon=False)
plt.show()
plt.close()

# %%
# https://github.com/tanghaibao/goatools/blob/main/notebooks/report_depth_level.ipynb
dict_GO_terms_up_hier_filt = dict()
for go_id, go_term in dict_GO_ID_terms_up.items():
    depth = go_depth.get(go_id, np.nan)
    if depth == 3:
        dict_GO_terms_up_hier_filt[go_term] = depth

dict_GO_terms_up_hier_filt

dict_GO_terms_down_hier_filt = dict()
for go_id, go_term in dict_GO_ID_terms_down.items():
    depth = go_depth.get(go_id, np.nan)
    if depth == 3:
        dict_GO_terms_down_hier_filt[go_term] = depth

dict_GO_terms_down_hier_filt

# %%
dict_up_bp_filt = {k: v for k,v in dict_up_bp.items() if k in dict_GO_terms_up_hier_filt.keys()}
dict_up_down_filt = {k: v for k,v in dict_down_bp.items() if k in dict_GO_terms_down_hier_filt.keys()}

df_up_top_filt = prepare_df(dict_up_bp_filt, "BP", top=10)
df_up_all_filt = prepare_df(dict_up_bp_filt, "BP", top=None)
df_down_top_filt = prepare_df(dict_up_down_filt, "BP", top=10)
df_down_all_filt = prepare_df(dict_up_down_filt, "BP", top=None)


from matplotlib.ticker import MultipleLocator

plt.rcParams["font.size"] = 14
fig, axs = plt.subplots(2, 1, figsize=(16, 10), facecolor="white")

sc_up = plot_lollipop(axs[0], df_up_top_filt, "Up-regulated", color="#d73027")
sc_down = plot_lollipop(axs[1], df_down_top_filt, "Down-regulated", color="#4575b4")

fig.text(0.02, 0.99, "A", fontsize=22, fontweight='bold', va='top', ha='left')
fig.text(0.02, 0.49, "B", fontsize=22, fontweight='bold', va='top', ha='left')

# Custom legend
legend_sizes = [1, 5, 10]
legend_dots = [
    plt.scatter([], [], s=s*20, edgecolors='black', facecolor='black', label=f"{s}", lw=0.5)
    for s in legend_sizes
]

# Add legends to each subplot
axs[0].legend(handles=legend_dots, title='-log10(qVal.)', loc='center', bbox_to_anchor=(1.5, 0.5), frameon=False)
axs[1].legend(handles=legend_dots, title='-log10(qVal.)', loc='center', bbox_to_anchor=(1.5, 0.5), frameon=False)

# Set x-axis tick spacing to 0.1
for ax in axs:
    ax.xaxis.set_major_locator(MultipleLocator(0.1))

axs[0].set_xlim(0, 0.4)
axs[1].set_xlim(0, 0.4)

plt.tight_layout()
plt.subplots_adjust(hspace=0.3)
# plt.savefig("")
plt.show()
plt.close()
# %%
