# %%
import subprocess

# %%
path_bbduk = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/bbmap/bbduk.sh"
in_fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.fasta"
out_fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.fasta"
adapter_fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/MGI_adapter_pe.fasta"

cmd = f"{path_bbduk} \
    in={in_fasta} \
    out={out_fasta} \
    ref={adapter_fasta} \
    ktrim=r \
    k=21 \
    rcomp=t \
    mink=11 \
    hdist=1 \
    tbo=t \
    tpe=t \
    minlen=20"

# %%
subprocess.run(cmd, shell=True)
