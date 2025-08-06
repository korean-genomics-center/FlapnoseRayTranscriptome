# %%
import os
import subprocess

## R requirements ##
# goseq (https://bioconductor.org/packages/release/bioc/html/goseq.html)
# qvalue (https://www.bioconductor.org/packages/release/bioc/html/qvalue.html)
# '''https://github.com/trinityrnaseq/trinity_ext_sample_data/tree/343f2eb856d3bcefbc84b4247b8109564ce9eae6/test_GOSeq_trinotate_pipe/Spombe_runGOseqOnly'''

# %%
def run_GOSeq(script, factor_label, go_annot, gene_lengths, expr_mtx):
    cmd = f"{script} --factor_labeling {factor_label} --GO_assignments {go_annot} --lengths {gene_lengths} --background {expr_mtx}"
    subprocess.run(cmd, shell=True)

# %%
script = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/trinityrnaseq-v2.14.0/Analysis/DifferentialExpression/run_GOseq.pl"
factor_label = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/GOSeq/factor_labeling_downreg_hit_cov_abov_60_20250725.tsv"
go_annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/GOSeq/trinotate_GO_annotation.txt"
gene_lengths = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.gene_lengths.txt"
expr_mtx = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.gene.TMM.EXPR.matrix.filt"

# %%
run_GOSeq(script, factor_label, go_annot, gene_lengths, expr_mtx)
