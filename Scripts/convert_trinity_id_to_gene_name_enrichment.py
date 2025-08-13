# %%
import os
from pathlib import Path


# %%
def get_hash_trinity_gene(annot):
    dict_trinity_gene = dict()
    with open(annot, "r") as fr:
        _skip = fr.readline()
        for line in fr:
            record = line.rstrip("\n").split("\t")
            trinity = record[0]
            gene_desc = record[2]
            if gene_desc == "NA":
                gene = "NA"
            else:
                gene = gene_desc.split(",")[1]
            dict_trinity_gene[trinity] = gene
    
    return dict_trinity_gene

def read_enrichment_result(enrich, hash, out):
    with open(out, "w") as fw:
        with open(enrich, "r") as fr:
            header = fr.readline()
            fw.write(header)
            for line in fr:
                list_gene_info = list()
                list_record = line.rstrip("\n").split("\t")
                list_trinity = list_record[-1].split(",")
                for trinity in list_trinity:
                    trinity = trinity.lstrip()
                    gene = hash[trinity]
                    gene_info = f"{trinity}:{gene}"
                    list_gene_info.append(gene_info)
                gene_info = ",".join(list_gene_info)
                list_record[-1] = gene_info
                record = "\t".join(list_record) + '\n'
                fw.write(record)

# %%
WORKDIR = str(Path(os.path.abspath(os.getcwd())).parents[1])
reg = "upreg"
annot = f"{WORKDIR}/Results/5_trinotate/trinotate_annotation_report_UniVec_fcs_filtered_att_hit_cov_abov_60.tsv"
enrich = f"{WORKDIR}/Results/7_go/diff_expressed_tumor_normal_{reg}.GOseq.enriched"
out = f"{WORKDIR}/Results/7_go/diff_expressed_tumor_normal_{reg}.GOseq.enriched_att_genename.tsv"

hash = get_hash_trinity_gene(annot)
read_enrichment_result(enrich, hash, out)
