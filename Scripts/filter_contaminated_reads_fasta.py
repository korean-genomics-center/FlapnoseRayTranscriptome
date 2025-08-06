# %%
import os
import pandas as pd
from fasta_parser import read_fasta_as_dict, save_dict_fasta_into_file
# %%
# path_UniVec = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/UniVec"
# dict_UniVec = read_fasta_as_dict(path_UniVec)
# dict_MGI_adapter = dict()
# for header, seq in dict_UniVec.items():
#     if "mgi" in header.lower():
#         dict_MGI_adapter[header] = seq

# path_MGI_adapter = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Resources/Database/MGI_adapter_pe.fasta"
# save_dict_fasta_into_file(dict_MGI_adapter, path_MGI_adapter)

# %%
path_fasta = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.fasta"
dict_fasta = read_fasta_as_dict(path_fasta)

# %%
# path_univec_blastn = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.outfmt6"
path_univec_blastn = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.second_run.outfmt6"
df_univect_blastn = pd.read_csv(path_univec_blastn, sep="\t", header=None)
set_univec_hit_trinity_id = set(df_univect_blastn[0].to_list())
set_univec_hit_trinity_id

# %%
dict_fasta_filtered_in = dict()
list_trinity_isoforms_filtered_in = list()
for header, seq in dict_fasta.items():
    trinity_id = header.split()[0]
    if trinity_id not in set_univec_hit_trinity_id and len(seq) > 200:
        list_trinity_isoforms_filtered_in.append(trinity_id)
        dict_fasta_filtered_in[header] = seq
        
# %%
dict_fasta_filtered_in

# %%
list_trinity_isoforms_filtered_in

# %%
list_trinity_genes_filtered_in = list(set(list(map(lambda x: "_".join(x.split("_")[:-1]), list_trinity_isoforms_filtered_in))))
list_trinity_genes_filtered_in

# %%
path_save = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Results/2_trinity/trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.fasta"
save_dict_fasta_into_file(dict_fasta_filtered_in, path_save)
# %%
