library(NOISeq)

path_exp <- "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/flapnoseray.tumor.normal.RSEM.genes.TMM.uniprot.annotated.hit.eval.1e05.coverage.abov.70.gct"
expr_ori <- read.csv(path_exp, sep="\t")
expr <- expr_ori[2:3]
rownames(expr_ori) <- as.matrix(expr_ori[1])
rownames(expr) <- rownames(expr_ori)

# Metadata
condition <- data.frame(condition = c("Left_tumor_inner", "normal"))
rownames(condition) <- colnames(expr)

# Prepare NOISeq object
mydata <- NOISeq::readData(data = expr, factors = condition)

# Run NOISeq with simulated technical replicates (N = 1 per group)
mynoiseq <- NOISeq::noiseq(
  input = mydata,
  factor = "condition",
  conditions = c("Left_tumor_inner", "normal"),
  replicates = "no",  # because you only have 1 sample per group
  nss = 3,            # number of simulated samples per condition
  norm = "n"          # no additional normalization (use "tmm" if you want)
)

# View results
results <- NOISeq::degenes(mynoiseq, q = 0.8)  # q = probability threshold

# Convert results to a data.frame
results_df <- as.data.frame(results)

# Calculate p-value and -log10(p-value)
results_df$p_value <- 1 - results_df$prob
results_df$neg_log10_pval <- -log10(results_df$p_value)

# Add FDR correction (Benjamini-Hochberg)
results_df$fdr <- p.adjust(results_df$p_value, method = "BH")

# Add protein labels
results_df$prot_label <- expr_ori[rownames(results_df), "gene_label"]

# View final table with combined label and FDR
results_df[, c("M", "D", "prob", "p_value", "fdr", "neg_log10_pval", "prot_label")]


path_deg <- "/BiO/Research/Project1/FlapnoserayTumorTranscriptome/Workspace/Resources/Scripts/FlapnoseRayTranscriptome/Data/NOISeq.deg.tsv"

write.table(results_df, path_deg, sep="\t", col.names=T, row.names=T, quote=F)