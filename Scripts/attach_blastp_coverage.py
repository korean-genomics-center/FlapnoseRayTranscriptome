#%%
from collections import Counter

# %%
def get_query_isoform_names(file):
    dict_query_cov = dict()
    with open(file, mode="r") as fr:
        for line in fr:
            record = line.rstrip("\n").split("\t")
            query_name = record[0].split(".")[0]
            hit_covgth = record[-2]
            dict_query_cov[query_name] = hit_covgth
    
    return dict_query_cov

def attach_blastp_cov(file, write_file, dict_query_cov, cov_thres):
    with open(write_file, "w") as fw:
        with open(file, mode="r") as fr:
            header = fr.readline().rstrip("\n").split("\t")
            new_header = "\t".join(header + ["hit_coverage"]) + "\n"
            fw.write(new_header)
            for line in fr: 
                record = line.rstrip("\n").split("\t")
                isoform = record[1]
                for query, hit_cov in dict_query_cov.items():
                    if isoform == query and float(hit_cov) > cov_thres:
                        line_att_hit_cov = "\t".join(record + [hit_cov]) + "\n"
                        fw.write(line_att_hit_cov)
                    else:
                        continue

def count_gene_isoform_number(path_new_report):
    list_gene = list()
    list_isoform = list()
    with open(path_new_report, mode="r") as fr:
        for line in fr:
            if line.startswith("TRINITY"):
                record = line.rstrip("\n").split("\t")
                gene = record[0]
                isoform = record[1]
                list_gene.append(gene)
                list_isoform.append(isoform)
            else:
                continue
    
    num_gene = len(dict(Counter(list_gene)))
    num_isoform = len(dict(Counter(list_isoform)))

    print("Number of genes:", num_gene)
    print("Number of isoform:", num_isoform)

# %%
cov_thres = 60
annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report_20250725.tsv"
cov = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/fulllengthtranscripts.w_pct_hit_length"
out = f"/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report_att_hit_cov_abov_{cov_thres}_20250725.tsv"

dict_query_cov = get_query_isoform_names(cov)
attach_blastp_cov(annot, out, dict_query_cov, cov_thres)
count_gene_isoform_number(out)

# %%
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(cov, sep="\t")
plt.figure(figsize=(5, 5), facecolor="white")
plt.axvline(x=cov_thres, linestyle="--", color="firebrick")
plt.rcParams["font.size"] = 14
plt.hist(df["pct_hit_len_aligned"], color="grey")
plt.xlabel("hit coverage (%)", fontsize=16)
plt.ylabel("gene count", fontsize=16)
plt.show()
plt.close()
# %%
