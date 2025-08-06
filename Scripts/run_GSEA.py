# %%
import pandas as pd
import gseapy as gp
from collections import defaultdict

# %%
def get_gmt_for_prerank(file):
    df = pd.read_csv(file, sep="\t", names=["gene", "GO_terms"])
    dict_GO_list_gene = defaultdict(list)
    for gene, GO_annot in zip(df["gene"], df["GO_terms"]):
        GO_terms = GO_annot.split(",")
        for GO in GO_terms:
            dict_GO_list_gene[GO].append(gene)   
    
    return dict_GO_list_gene 

# %%
file = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/6_deg/trinotate_GO_annotation.txt"
gmt = get_gmt_for_prerank(file)

# %%
path_log2fc = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.rnk"
ranking_df = pd.read_csv(path_log2fc, sep="\t")

# %%
# Step 1: Prepare a ranked list — a dataframe with 'gene_name' and 'ranking_score'
ranking_df = ranking_df[['index', 'log2fc']]
ranking_df.columns = ['gene_name', 'log2fc']
ranking_list = ranking_df.set_index("gene_name")["log2fc"].sort_values(ascending=False)

# %%
# # Step 2: Run GSEA using GO terms (downloaded or prebuilt)
# gsea_results = gp.prerank(
#     rnk=ranking_list,
#     gene_sets="gmt",
#     processes=4,
#     permutation_num=1000,
#     outdir='/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results',
#     seed=6,
#     format='png',
# )

# %%
# Step 2: Run GSEA using GO terms (downloaded or prebuilt)
gsea_results = gp.prerank(
    rnk=ranking_list,
    gene_sets=gmt,
    processes=4,
    permutation_num=1000,
    outdir='/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results',
    seed=6,
    format='png',
)

# %%
# Step 3: View results
gsea_sig = gsea_results.res2d[gsea_results.res2d["FDR q-val"] < 0.05]

# %%
