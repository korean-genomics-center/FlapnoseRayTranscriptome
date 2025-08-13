# %%
import os
from pathlib import Path

import pandas as pd

from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])

# %%
path_fasta = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.contam_cleaned.fasta"
dict_fasta = read_fasta_as_dict(path_fasta)
trinity_isoforms = list(map(lambda x: x.split()[0].strip(), list(dict_fasta.keys())))

# %%
path_genetransmap = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map"
dict_transgenemap = dict()
with open(path_genetransmap, mode="r") as fr:
    for line in fr:
        record = line.rstrip().split()
        gene = record[0]
        trans = record[1]
        dict_transgenemap[trans] = gene 

trinity_genes = list(set(list(map(lambda x: dict_transgenemap.get(x, None), trinity_isoforms))))

# %%
path_annot = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report.xls"
df_annot = pd.read_csv(path_annot, sep="\t")
df_annot_select = df_annot[df_annot["transcript_id"].isin(trinity_isoforms)]

path_annot_filt = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered.xls"
df_annot_select.to_csv(path_annot_filt, sep="\t", index=False)

# %%
exp_matrix = f"{WORKDIR}/Results/3_rsem/RSEM.isoform.TMM.EXPR.matrix"
df_exp = pd.read_csv(exp_matrix, sep="\t").rename(columns={"Unnamed: 0": "transcript_id"})
df_exp_select = df_exp[df_exp["transcript_id"].isin(trinity_isoforms)]
df_exp_select = df_exp_select.set_index("transcript_id")
df_exp_select.index.name = None

exp_matrix_filt = f"{WORKDIR}/Results/3_rsem/RSEM.isoform.TMM.EXPR.matrix.filt"
df_exp_select.to_csv(exp_matrix_filt, sep="\t", index=True)


# %%
exp_matrix = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix"
df_exp = pd.read_csv(exp_matrix, sep="\t").rename(columns={"Unnamed: 0": "gene_id"})
df_exp_select = df_exp[df_exp["gene_id"].isin(trinity_genes)]
df_exp_select = df_exp_select.set_index("gene_id")
df_exp_select.index.name = None

exp_matrix_filt = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"
df_exp_select.to_csv(exp_matrix_filt, sep="\t", index=True)

# %%
gene_trans_map = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map"
df_gene_gene_trans_map = pd.read_csv(gene_trans_map, sep="\t", header=None, names=["gene_id", "transcript_id"])
df_gene_gene_trans_map_select = df_gene_gene_trans_map[df_gene_gene_trans_map["transcript_id"].isin(trinity_isoforms)]

df_gene_gene_trans_map_filt = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map.filt"
df_gene_gene_trans_map_select.to_csv(df_gene_gene_trans_map_filt, sep="\t", index=False, header=False)
