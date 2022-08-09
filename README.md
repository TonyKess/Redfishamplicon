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
ls * ls *R1_001.fastq.gz | \
  sed 's/_R1_001.fastq.gz//' | 
  sort |
  uniq > inds.tsv
```

And align the reads

```
while read ind;
  do echo $ind\.bam ;
  RGID=$(echo $ind |  sed 's/i5.*/i5/') ;
  SMID=$(echo $ind | sed 's/NS.*i5.//') ;
  LBID=$(echo $ind | sed 's/.UDP.*//');
  $bwamem2 mem \
  -t 32 \
  -R "@RG\tID:$RGID\tSM:$SMID\tLB:$LBID" \
  $genome \
  $ind\_R1.trimmed.fastq.gz  $ind\_R2.trimmed.fastq.gz\
  | samtools sort -o $projdir/align/$ind\.sorted.bam -T $ind -@ 32 -m 3G ;
  done <  inds.tsv
```
