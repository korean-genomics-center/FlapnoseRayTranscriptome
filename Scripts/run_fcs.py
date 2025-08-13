# %%
import os
import subprocess

# %%
fcs = "/BiO/Access/kyungwhan1998/FlapnoserayTumorTranscriptome/Resources/Tools/fcs-0.5.5/dist/fcs.py"
path_fa = "/BiO/Access/kyungwhan1998/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.fasta"
gx_out = "/BiO/Access/kyungwhan1998/FlapnoserayTumorTranscriptome/Results/8_fcs"
gx_db = "/BiO/Access/kyungwhan1998/FlapnoserayTumorTranscriptome/Resources/Database/gxdb"
tax_id = "195332"

cmd = f"{fcs} screen genome --fasta {path_fa} --out-dir {gx_out} --gx-db {gx_db} --tax-id {tax_id}"

# %%
subprocess.run(cmd, shell=True)