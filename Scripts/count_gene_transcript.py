#%%

from collections import Counter

annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/6_deg/trinotate_annotation_report_20230619.tsv"
list_gene = list()
list_isoform = list()
with open(annot, 'r') as fr:
    _skiprow = fr.readline()
    for line in fr:
        record = line.rstrip("\n").split("\t")
        gene = record[0]
        isoform = record[1]
        list_gene.append(gene)
        list_isoform.append(isoform)

count_gene = len(dict(Counter(list_gene)).keys())
count_isoform = len(dict(Counter(list_isoform)).keys())

print(count_gene)
print(count_isoform)

#%%

from collections import Counter

annot = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/5_trinotate/trinotate_annotation_report.xls"
list_gene = list()
list_isoform = list()
with open(annot, 'r') as fr:
    _skiprow = fr.readline()
    for line in fr:
        record = line.rstrip("\n").split("\t")
        pep = record[4]
        if pep == ".":
            continue
        gene = "_".join(pep.split("_")[:4])
        isoform = "_".join(pep.split("_")[:5])
        list_gene.append(gene)
        list_isoform.append(isoform)

count_gene = len(dict(Counter(list_gene)).keys())
count_isoform = len(dict(Counter(list_isoform)).keys())

print(count_gene)
print(count_isoform)