# FlapnoseRayTranscriptome

### Content
* [Analysis](#analysis)
  - [Pre-Assembly Processing and Quality Control](#pre-assembly-processing-and-quality-control)
  - [De novo Transcriptome Assembly](#merge-cpg-values-into-single-table)
  - [Post-assembly Quality Control](#run-weighted-differential-methylation-analysis)
  - [ORF Prediction](#orf-prediction)
  - [Differential Expression Analysis](#differential-expression-analysis)
  - [Functional Enrichment](#functional-enrichment)

## Analysis
### Pre-Assembly Processing and Quality Control
* Read trimming ([fastp](https://github.com/OpenGene/fastp))

```bash
fastp --thread 16 --in1 Left_tumor_inner_1.fq.gz --in2 Left_tumor_inner_2.fq.gz --out1 Left_tumor_inner_1.trimmed.fq.gz --out2 Left_tumor_inner_2.trimmed.fq.gz --json Left_tumor_inner.json --html Left_tumor_inner.html --report_title Left_tumor_inner --cut_front 1 --cut_right 0 --cut_tail 1 --detect_adapter_for_pe --trim_front1 15 --trim_front2 15 --trim_tail1 5 --trim_tail2 5 --cut_mean_quality 20 --n_base_limit 1 --average_qual 20
```

### De novo Transcriptome Assembly
* Transcriptome Assembly ([Trinity](https://github.com/trinityrnaseq/trinityrnaseq))

```bash
Trinity --seqType fq --left Left_tumor_1.fq.gz, Normal_1.fq.gz --right Left_tumor_2.fq.gz Normal_2.fq.gz --CPU 70 –max_memory 100G --output trinity_flapnoseray
```

### Post-assembly Quality Control
* Check BUSCO completeness ([BUSCO](https://busco.ezlab.org/busco_userguide.html))

```bash
busco -c 30 -i trinity_flapnoseray.Trinity.fasta --out_path busco -o {lineage} -m transcriptome -f 
```

* Remove remaining adapters from read ([bbduck](https://sourceforge.net/projects/bbmap/))

```bash
bbduk.sh \
    in=trinity_flapnoseray.Trinity.fasta \
    out=trinity_flapnoseray.Trinity.bbuk.adapter.filtered.fasta \
    ref=MGI_adapter_pe.fasta \
    ktrim=r \
    k=21 \
    rcomp=t \
    mink=11 \
    hdist=1 \
    tbo=t \
    tpe=t \
    minlen=20
```

* Screen UniVec database ([UniVec](https://www.ncbi.nlm.nih.gov/genbank/tsafaq/))

```bash
makeblastdb -in UniVec -dbtype nucl -out UniVec

blastn -task blastn \
    -reward 1 \
    -penalty -5 \
    -gapopen 3 \
    -gapextend 3 \
    -dust yes \
    -soft_masking true \
    -evalue 700 \
    -searchsp 1750000000000 \
    -db UniVec \
    -query trinity_flapnoseray.Trinity.bbuk.adapter.filtered.fasta \
    -outfmt 6 \
    -out trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.outfmt6
```

* Foreign contamination screening ([FCS-GX](https://github.com/ncbi/fcs-gx))

```bash
python fcs.py screen genome --fasta trinity_flapnoseray.Trinity.bbuk.adapter.filtered.UniVec.blastn.filtered.fasta --out-dir gx_out --gx-db gx_db --tax-id 195332
```

* Remove spurius genes/transcripts 
```bash
python remove_spurious_genes.py
```

### ORF Prediction

* TransDecoder LongOrfs ([TransDecoder](https://github.com/TransDecoder/TransDecoder))

```bash
TransDecoder.LongOrfs -t trinity_flapnoseray.Trinity.fasta --gene_trans_map trinity_flapnoseray.Trinity.fasta.gene_trans_map -m 100 --output_dir transdecoder --complete_orfs_only
```

* TransDecoder Predict ([TransDecoder](https://github.com/TransDecoder/TransDecoder))

```bash
TransDecoder.Predict -t trinity_flapnoseray.Trinity.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only 
```

* Analyze full-length transcript ([analyze_blastPlus_topHit_coverage.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Counting-Full-Length-Trinity-Transcripts))

```bash
analyze_blastPlus_topHit_coverage.pl blastp.outfmt6 trinity_flapnoseray.Trinity.fasta uniprot_sprot.fasta
```

* Eliminate redundant transcripts ([CD-HIT](https://bioinformatics.org/cd-hit/))

#### first pass

```bash
cd-hit-est -T 30 -c 0.98 -M 0 -d 0 -b 3 -p 1 \
-i trinity_flapnoseray.Trinity.fasta.transdecoder.pep \                      
-o trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass1_cdhit0.98
```
#### second pass
```bash
cd-hit-est -T 30 -c 0.98 -M 0 -d 0 -b 3 -p 1 \
-i trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass1_cdhit0.98 \
-o trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass2_cdhit0.98
```

* Annotate Protein-coding Transcripts ([Trinotate](https://github.com/Trinotate/Trinotate))

```bash 
Trinotate --db flapnoseray.sqlite init --gene_trans_map trinity_flapnoseray.Trinity.fasta.gene_trans_map --transcript_fasta trinity_flapnoseray.Trinity.fasta --transdecoder_pep trinity_flapnoseray.Trinity.fasta.transdecoder.pep_pass2_cdhit0.98. 

Trinotate --db flapnoseray.sqlite --LOAD_swissprot_blastp blastp.outfmt6 --LOAD_doublout_pfam pfam.domtblout

Trinotate --db flapnoseray.sqlite --report -E 1e-05 > trinotate_annotation_report.xls
```

## Differential Expression Analysis 

* Build bowtie2 reference ([Bowtie2](https://github.com/Trinotate/Trinotate))
```bash
bowtie2-build –threads 30 trinity_flapnoseray.Trinity.fasta trinity_flapnoseray.Trinity.fasta.bowtie2
```

* Build RSEM reference ([RSEM](https://github.com/Trinotate/Trinotate))
```bash
rsem-prepare-reference -p 30 trinity_flapnoseray.Trinity.fasta --transcript-to-gene-map trinity_flapnoseray.Trinity.fasta.gene_trans_map --bowtie2 --bowtie2-path bowtie2 trinity_flapnoseray.Trinity.fasta.RSEM
```

* Align reads to transcriptome assembly ([Bowtie2](https://github.com/Trinotate/Trinotate))
```bash
bowtie2 -p 50 --no-unal --no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 -X 800 -x trinity_flapnoseray.Trinity.fasta.RSEM -q -1 Left_tumor_1.fq.gz (or Normal_1.fq.gz) -2 Left_tumor_2.fq.gz (or Normal_2.fq.gz ) 2> align_stats.txt
```

* Quantify read counts ([RSEM](https://github.com/Trinotate/Trinotate))
```bash
rsem-calculate-expression -p 8 --paired-end –alignments --estimate-rspd --append-names --no-bam-output trinity_flapnoseray.Trinity.fasta.bowtie2.bam trinity_flapnoseray.Trinity.fasta.RSEM  Left_tumor(or Normal)/RSEM
```

* Create expression matrix ([abundance_estimates_to_matrix.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification))
```bash
abundance_estimates_to_matrix.pl --out_prefix RSEM --gene_trans_map trinity_flapnoseray.Trinity.fasta.gene_trans_map --cross_sample_norm TMM --est_method RSEM --name_sample_by_basedir Left_tumor_inner/RSEM.isoforms.results normal/RSEM.isoforms.results
```

## Functional Enrichment
* Make factor label
```bash
python make_factor_labeling.py
```

* Get gene length ([fasta_seq_length.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq))
```bash
python get_tpm_weighted_gene_length_trinity.py
```

* Extract GO assignment ([extract_GO_assignments_from_Trinotate_xls.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq))

```bash
python extract_GO_assignment.py
```

* Gene ontology enrichments ([runGOSeq.pl](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-GOSeq))

```bash
python run_GOSeq.py
```
