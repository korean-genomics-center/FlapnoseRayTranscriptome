# %%
import subprocess

TRINOTATE = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Tools/Trinotate-Trinotate-v4.0.1"
annotation_file = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report.xls"
GO_file = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/GOSeq/trinotate_GO_annotation.txt"

cmd = f"{TRINOTATE}/util/extract_GO_assignments_from_Trinotate_xls.pl \
                         --Trinotate_xls {annotation_file} \
                         -G --include_ancestral_terms \
                         > {GO_file}"

subprocess.run(cmd, shell=True)
