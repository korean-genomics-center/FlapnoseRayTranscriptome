# %%
def get_cnt_lines(file):
    cnt = 0 
    with open(file, mode="r") as fr:
        for line in fr:
            if str(line).startswith(">"):
                cnt += 1

    return cnt

# %%
long = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/4_transdecoder/longest_orfs.pep"
cnt_long = get_cnt_lines(long)
pred = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/4_transdecoder/trinity_flapnoseray.Trinity.fasta.transdecoder.pep"
cnt_pred = get_cnt_lines(pred)
pass1 = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/4_transdecoder/trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass1_cdhit0.98"
cnt_pass1 = get_cnt_lines(pass1)
pass2 = "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Results/20230701/4_transdecoder/trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass2_cdhit0.98"
cnt_pass2 = get_cnt_lines(pass2)

print(cnt_long)
print(cnt_pred)
print(cnt_pass1)
print(cnt_pass2)
