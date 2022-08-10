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

Index the genome

```
bwa index genome/GCA_916700875.1_S-aleutianus_SEB-111_genomic.fna
```

Make a list of reads

```
cd amplicon_fastq

ls trim/*.fastq.gz  | \
  sed 's/trim\///'  | \
  sed 's/.trim.*//' | \
  sort | uniq > inds.tsv
```

And align the reads

```
mkdir align

while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind) ;
  SMID=$(echo $ind) ;
  LBID=$(echo $ind) ;
  bwa mem \
  -t 32 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
   ../genome/GCA_916700875.1_S-aleutianus_SEB-111_genomic.fna \
  trim/$ind.trimmed.fastq.gz | samtools sort -o align/$ind\.sorted.bam -T $ind -@ 24  ;
done <  inds.tsv
```

Deduplicate reads

```
#deduplicate in parallel
cat inds.tsv | \
  parallel --jobs 24 'gatk --java-options "-Xmx4G" \
  MarkDuplicates \
  I=align/{}.sorted.bam \
  O=align/{}.deDup.bam M=align/{}_deDupMetrics.txt \
  REMOVE_DUPLICATES=true'
```

Index after deduplication
```
while read  ind; 
  do samtools index align/$ind.deDup.bam ; done < inds.tsv
```

Make sequence dictionary
```
gatk CreateSequenceDictionary -R ../genome/GCA_916700875.1_S-aleutianus_SEB-111_genomic.fna
```


Haplotype caller genotyping
```
mkdir vcfs 

ls align/*.deDup.bam | sed 's/.deDup.bam//' | \
  parallel --jobs 24 \
  'gatk --java-options "-Xmx2G" \
  HaplotypeCaller \
  -ERC GVCF \
  -R ../genome/GCA_916700875.1_S-aleutianus_SEB-111_genomic.fna \
  -I align/{}.deDup.bam -O vcfs/{}.g.vcf'
```
