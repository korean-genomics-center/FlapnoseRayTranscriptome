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
reg = "upreg"
annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Data/Trinotate/trinotate_annotation_report_att_hit_cov_abov_60_20250725.tsv"
enrich = f"/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Scripts/diff_expressed_tumor_normal_{reg}.GOseq.enriched"
out = f"/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/FlapnoseRayTranscriptome/Scripts/diff_expressed_tumor_normal_{reg}.GOseq.enriched_att_genename.tsv"

hash = get_hash_trinity_gene(annot)
read_enrichment_result(enrich, hash, out)