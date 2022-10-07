# Redfishamplicon

Make a conda environment for fastp, activate, trim reads in parallel

```
conda create -n fastp fastp -c bioconda

conda activate fastp

ls *R1_001.fastq.gz | \
  sed 's/_R1_001.fastq.gz//' | parallel --jobs 8 \
  'fastp -i {}_R1_001.fastq.gz \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  -o {}.trimmed.fastq.gz \
  --thread 4 '
  
  conda deactivate
  ```

Index the amplicon

```
bwa index fasciatus_mentella_LE_fst_top50.fasta
```

Align to amplicon

```
ls trim/*trimmed.fastq.gz | \
  sed 's/.trimmed.fastq.gz//' | \
  sed 's/trim\///' | \
  parallel \
  -j 8 \
    'bwa mem -R "@RG\tID:{}\tSM:{}\tLB:library1" \
    fasciatus_mentella_LE_fst_top50.fasta trim/{}.trimmed.fastq.gz | \
    samtools sort -o amplicon_align/{}.sorted.bam '
```
SNP call in bcftools

```
bcftools mpileup \
  -Ou \
  --max-depth 2000000 \
  --threads 24  \
  -f fasciatus_mentella_LE_fst_top50.fasta \
 amplicon_align/*sorted.bam | \
    bcftools call \
    --threads 24 \
    -mv \
    -Ob \
    -V indels \
    -o redfishamplicon.bcf
```    
Filter by depth (180 x 10 and qual 30)

```
 
bcftools filter -i "DP>1800  && QUAL >= 30" redfishamplicon.bcf  > redfishamplicon_dp1K_q30.vcf

```
Filter for missing and remove unknowns in plink

```
plink --vcf redfishamplicon_dp1K_q30.vcf --allow-extra-chr --make-bed --double-id --maf 0.01 --geno 0.05 --mind 0.05 --out redfishamplicon_dp1K_q30_geno05 --remove no_UNK
```
Get some basic IDs for plotting
```
awk '{print $1}' redfishamplicon_dp1K_q30_geno05.fam | sed 's/\-.*//' > amp_IDs
```
