# %%
import subprocess
import os

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
get_gene_length = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/trinityrnaseq-v2.14.0/util/misc/fasta_seq_length.pl"
get_weighted_gene_length = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/trinityrnaseq-v2.14.0/util/misc/TPM_weighted_gene_length.py"
fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.fasta"
gene_trans_map = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.fasta.gene_trans_map"
exp_matrix = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/3_rsem/RSEM.isoform.TMM.EXPR.matrix"

# %%
seq_lens, cmd1 = make_gene_length_cmd(get_gene_length, fasta)
subprocess.run(cmd1, shell=True)

# %%
cmd2 = make_weighted_gene_length_cmd(get_weighted_gene_length, gene_trans_map, seq_lens, exp_matrix)
subprocess.run(cmd2, shell=True)
# %%
