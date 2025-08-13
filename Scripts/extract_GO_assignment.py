# %%
import os
import subprocess
from pathlib import Path

# %%

WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
TRINOTATE = f"{WORKDIR}/Resources/Tools/Trinotate-Trinotate-v4.0.1"
annotation_file = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered.xls"
GO_file = f"{WORKDIR}/Results/7_go/trinotate_GO_annotation_UniVec_fcs_filtered.txt"

cmd = f"{TRINOTATE}/util/extract_GO_assignments_from_Trinotate_xls.pl \
        --Trinotate_xls {annotation_file} \
        -G --include_ancestral_terms \
        > {GO_file}"

subprocess.run(cmd, shell=True)
