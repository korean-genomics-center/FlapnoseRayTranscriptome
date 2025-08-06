# %%
import os
import pandas as pd
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file

# %%
path_fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.fasta"
dict_fasta = read_fasta_as_dict(path_fasta)
trinity_isoforms = list(map(lambda x: x.split()[0].strip(), list(dict_fasta.keys())))

# %%oforms

path_genetransmap = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map"
dict_transgenemap = dict()
with open(path_genetransmap, mode="r") as fr:
    for line in fr:
        record = line.rstrip().split()
        gene = record[0]
        trans = record[1]
        dict_transgenemap[trans] = gene 

trinity_genes = list(set(list(map(lambda x: dict_transgenemap.get(x, None), trinity_isoforms))))

# %%
path_annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report.xls"
df_annot = pd.read_csv(path_annot, sep="\t")
df_annot_select = df_annot[df_annot["transcript_id"].isin(trinity_isoforms)]

# path_annot_filt = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report_UniVec_filtered.xls"
# df_annot_select.to_csv(path_annot_filt, sep="\t", index=False)

# %%
exp_matrix = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix"
df_exp = pd.read_csv(exp_matrix, sep="\t").rename(columns={"Unnamed: 0": "gene_id"})
df_exp_select = df_exp[df_exp["gene_id"].isin(trinity_genes)]
df_exp_select = df_exp_select.set_index("gene_id")
df_exp_select.index.name = None

exp_matrix_filt = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"
df_exp_select.to_csv(exp_matrix_filt, sep="\t", index=True)