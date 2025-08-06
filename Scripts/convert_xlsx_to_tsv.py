import pandas as pd


path_original_report = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/6_deg/20230703/trinotate_annotation_report_20230703.tsv"
path_high_hit_cov_report = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/6_deg/20230703/trinotate_annotation_report_att_hit_cov_20230703.tsv"
path_enrich_upreg = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/6_deg/20230703/diff_expressed_tumor_normal_upreg.GOseq.enriched_att_genename.tsv"
path_enrich_downreg = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/6_deg/20230703/diff_expressed_tumor_normal_downreg.GOseq.enriched_att_genename.tsv"

def convert_tsv_xlsx(tsv):
    xlsx = tsv.split(".")[0] + ".xlsx"
    df = pd.read_csv(tsv, sep="\t")
    df.to_excel(xlsx, index=False)

# convert_tsv_xlsx(path_original_report)
# convert_tsv_xlsx(path_high_hit_cov_report)
convert_tsv_xlsx(path_enrich_upreg)
convert_tsv_xlsx(path_enrich_downreg)
