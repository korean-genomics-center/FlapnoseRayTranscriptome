# %%
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


# %%
def read_transcript_annot_report(annot_report):
    df_annot = pd.read_csv(annot_report, sep="\t")
    df_annot = df_annot[["#gene_id", "transcript_id", "sprot_Top_BLASTP_hit", "Pfam", "gene_ontology_BLASTP"]]
    return df_annot

def get_blastp_result(x):
    try:
        gene_sym = x.split("`")[0].split("^")[0].split("_")[0]
        organism = x.split("`")[0].split("^")[0].split("_")[1]
        gene_name = x.split("`")[0].split("^")[5].split("Full")[1][1:-1].split("{")[0]
        e_value = x.split("`")[0].split("^")[4]
        res = f"{gene_name},{gene_sym},{organism},{e_value}"
        return res
    except:
        return "NA"

def get_pfam_result(x):
    try:
        protein_name = x.split("^")[1]
        e_value = x.split("^")[-1]
        return f"{protein_name},{e_value}"
    except:
        return "NA"

def get_go_result(x):
    try:
        go = x.split("`")[0].split("^")[0]
        func = x.split("`")[0].split("^")[-1]
        return f"{go},{func}"
    except:
        return "NA"

def modify_annot_columns(df_annot):
    df_annot = df_annot[(df_annot["Pfam"]!=".")|(df_annot["sprot_Top_BLASTP_hit"]!=".")]
    df_annot["sprot_Top_BLASTP_hit"] = df_annot["sprot_Top_BLASTP_hit"].apply(get_blastp_result)
    df_annot["Pfam"] = df_annot["Pfam"].apply(get_pfam_result)
    df_annot["gene_ontology_BLASTP"] = df_annot["gene_ontology_BLASTP"].apply(get_go_result)
    df_annot.columns = ["GENE", "ISOFORM", "BLASTP", "PFAM", "GO"]
    return df_annot

def __read_exp_matrix(exp_matrix):
    df_exp = pd.read_csv(exp_matrix, sep="\t")
    df_exp["log2FC"] = np.log2((df_exp["Left_tumor_inner"] + 1) / (df_exp["normal"] + 1))
    df_exp.columns = ["gene_id", "tumor", "normal", "log2FC"]
    df_exp = df_exp.sort_values(by="log2FC", key=abs, ascending=False)
    return df_exp

def get_dict_query_exp(exp_matrix):
    df_exp = __read_exp_matrix(exp_matrix)
    return dict(zip(df_exp["gene_id"], df_exp["log2FC"].astype(str)))

def write_annot_file(df_annot, dict_query_exp, fileout):
    with open(fileout, "w") as fw:
        header = "\t".join(["GENE", "ISOFORM", "BLASTP", "PFAM", "GO", "LOG2FC"]) + "\n"
        fw.write(header)
        for gene, isoform, blastp, pfam, go in zip(df_annot["GENE"], df_annot["ISOFORM"], df_annot["BLASTP"], df_annot["PFAM"], df_annot["GO"]):
            if gene in dict_query_exp:
                log2fc = dict_query_exp[gene]
                content = "\t".join([gene, isoform, blastp, pfam, go, log2fc]) + "\n"
                fw.write(content)

def count_gene_isoform_number(path_new_report):
    list_gene = []
    list_isoform = []
    with open(path_new_report, mode="r") as fr:
        for line in fr:
            if line.startswith("TRINITY"):
                record = line.rstrip("\n").split("\t")
                gene = record[0]
                isoform = record[1]
                list_gene.append(gene)
                list_isoform.append(isoform)
    num_gene = len(dict(Counter(list_gene)))
    num_isoform = len(dict(Counter(list_isoform)))
    print("Number of genes:", num_gene)
    print("Number of isoforms:", num_isoform)

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
annot_report = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered.xls"
exp_matrix = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"
file_out = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered.tsv"
os.makedirs(os.path.dirname(file_out), exist_ok=True)

# %%
df_annot = read_transcript_annot_report(annot_report)
df_annot_modif = modify_annot_columns(df_annot)
dict_query_exp = get_dict_query_exp(exp_matrix)
write_annot_file(df_annot_modif, dict_query_exp, file_out)
count_gene_isoform_number(file_out)

# %%
