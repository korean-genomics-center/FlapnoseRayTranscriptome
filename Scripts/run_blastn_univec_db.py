# %%
import os
import subprocess

# %%
path_makeblastdb = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/ncbi-blast-2.13.0+/bin/makeblastdb"
path_univec_faa = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/UniVec"
makeblastdbout = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/UniVec"

cmd = f"{path_makeblastdb} -in {path_univec_faa} -dbtype nucl -out {makeblastdbout}"
# subprocess.run(cmd, shell=True)

 # %%
path_blastn = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/ncbi-blast-2.13.0+/bin/blastn"
path_trinity_faa = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.fasta"
blastnout = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.outfmt6"

blastn_cmd = f"{path_blastn} \
    -task blastn \
    -reward 1 \
    -penalty -5 \
    -gapopen 3 \
    -gapextend 3 \
    -dust yes \
    -soft_masking true \
    -evalue 700 \
    -searchsp 1750000000000 \
    -db {makeblastdbout} \
    -query {path_trinity_faa} \
    -outfmt 6 \
    -out {blastnout}"

# subprocess.run(blastn_cmd, shell=True)

# %%
path_blastn = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/ncbi-blast-2.13.0+/bin/blastn"
path_trinity_faa = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.fasta"
blastnout = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.second_run.outfmt6"

blastn_cmd = f"{path_blastn} \
    -task blastn \
    -reward 1 \
    -penalty -5 \
    -gapopen 3 \
    -gapextend 3 \
    -dust yes \
    -soft_masking true \
    -evalue 700 \
    -searchsp 1750000000000 \
    -db {makeblastdbout} \
    -query {path_trinity_faa} \
    -outfmt 6 \
    -out {blastnout}"

subprocess.run(blastn_cmd, shell=True)
# %%
