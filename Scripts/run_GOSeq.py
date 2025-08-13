# %%
import os
import subprocess
from pathlib import Path

## R requirements ##
# goseq (https://bioconductor.org/packages/release/bioc/html/goseq.html)
# qvalue (https://www.bioconductor.org/packages/release/bioc/html/qvalue.html)
# '''https://github.com/trinityrnaseq/trinity_ext_sample_data/tree/343f2eb856d3bcefbc84b4247b8109564ce9eae6/test_GOSeq_trinotate_pipe/Spombe_runGOseqOnly'''

# %%
def run_GOSeq(script, factor_label, go_annot, gene_lengths, expr_mtx):
    cmd = f"{script} --factor_labeling {factor_label} --GO_assignments {go_annot} --lengths {gene_lengths} --background {expr_mtx}"
    subprocess.run(cmd, shell=True)

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
script = f"{WORKDIR}/Resources/Tools/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_GOseq.pl"
factor_label = f"{WORKDIR}/Results/7_go/factor_labeling_upreg_UniVec_fcs_filtered_att_hit_cov_abov_60.tsv"
go_annot = f"{WORKDIR}/Results/7_go/trinotate_GO_annotation_UniVec_fcs_filtered.txt"
gene_lengths = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.contam_cleaned.gene_lengths.txt"
expr_mtx = f"{WORKDIR}/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt.lowcount.removed"

# %%
run_GOSeq(script, factor_label, go_annot, gene_lengths, expr_mtx)

# %%
