# %%
import os
import subprocess
from pathlib import Path


# %%
def make_gene_length_cmd(script, fasta):
    outfile = fasta.split(".fasta")[0] + ".seq_lens"
    cmd = f"{script} {fasta} > {outfile}"

    return outfile, cmd

def make_weighted_gene_length_cmd(script, gene_trans_map, seq_lens, exp_matrix):
    outfile = seq_lens.split(".seq_lens")[0] + ".gene_lengths.txt"
    cmd = f"{script} --gene_trans_map {gene_trans_map} --trans_lengths {seq_lens} --TPM_matrix {exp_matrix} > {outfile}"

    return cmd

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
get_gene_length = f"{WORKDIR}/Resources/Tools/trinityrnaseq-v2.14.0/util/misc/fasta_seq_length.pl"
get_weighted_gene_length = f"{WORKDIR}/Resources/Tools/trinityrnaseq-v2.14.0/util/misc/TPM_weighted_gene_length.py"
fasta = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.contam_cleaned.fasta"
gene_trans_map = f"{WORKDIR}/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map.filt"
exp_matrix = f"{WORKDIR}/Results/3_rsem/RSEM.isoform.TMM.EXPR.matrix.filt"

# %%
seq_lens, cmd1 = make_gene_length_cmd(get_gene_length, fasta)
subprocess.run(cmd1, shell=True)

# %%
cmd2 = make_weighted_gene_length_cmd(get_weighted_gene_length, gene_trans_map, seq_lens, exp_matrix)
subprocess.run(cmd2, shell=True)
# %%
