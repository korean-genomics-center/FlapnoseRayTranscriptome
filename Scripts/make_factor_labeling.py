# %%
import os
from pathlib import Path


# %%
def get_unique_deg(infile, reg):
    seen_genes = set()
    with open(infile, mode='r') as fr:
        _ = fr.readline()  # skip header
        for line in fr:
            record = line.rstrip().split("\t")
            if len(record) < 2:
                continue
            gene_id = record[0]
            log2fc = record[-2]
            try:
                log2fc_val = float(log2fc)
            except ValueError:
                continue

            if reg == "downreg" and log2fc_val < 0 and gene_id not in seen_genes:
                seen_genes.add(gene_id)
            elif reg == "upreg" and log2fc_val > 0 and gene_id not in seen_genes:
                seen_genes.add(gene_id)
    return list(seen_genes)

def get_sig_deg(outfile, list_deg, reg):
    with open(outfile, mode='w') as fw:
        for deg in list_deg:
            fw.write(f"diff_expressed_tumor_normal_{reg}\t{deg}\n")

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
reg = "downreg"
infile = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered_att_hit_cov_abov_60.tsv"
outfile = f"{WORKDIR}/Results/7_go/factor_labeling_{reg}_UniVec_fcs_filtered_att_hit_cov_abov_60.tsv"
list_deg = get_unique_deg(infile, reg)
get_sig_deg(outfile, list_deg, reg)
