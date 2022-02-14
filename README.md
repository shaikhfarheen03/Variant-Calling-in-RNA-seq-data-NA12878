# Variant Calling in NA12878 RNA data
Variant Calling in NA12878 was performed using the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-).
 
 # Fastq files and reference genome
Fastq files were obtained from the following link https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5665260 
To download the files I installed SRAtoolkit on the server and used the fasterq dump command: 
```bash
fasterq-dump SRR5665260
```

Reference genome was downloaded from the following link: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false 
The fasta file, .fai file, .dict, dbsnp file, known indel sites were also downloaded.

GTF file was downloaded from Gencode: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
 
 
 # Workflow
 ![Workflow](https://user-images.githubusercontent.com/26681884/153778754-a7d5fcd5-69f2-4ac5-8fa2-b0cc76198775.jpg) 
 [RNA variant calling steps and parameters - Sheet1.pdf](https://github.com/shaikhfarheen03/Varaint-Calling-in-RNA-seq-data-NA12878-/files/8056851/RNA.variant.calling.steps.and.parameters.-.Sheet1.pdf)

 ## Running fastp for quality check and to trim adapters 
Input 
>SRR5665260_1.fastq, \
>SRR5665260_2.fastq

Command
```bash 
fastp -i SRR5665260_1.fastq -I SRR5665260_2.fastq -o SRR5665260_1.fastp.fastq -O SRR5665260_2.fastp.fastq --detect_adapter_for_pe 
```
Output
>SRR5665260_1.fastp.fastq, \
>SRR5665260_2.fastp.fastq, \
>fastp.html, \
>fastp.json

 ## Generating a index for the reference genome using the gtf file
Input 
>resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta,# reference genome\
>gencode_v39_annotation.gtf #gtf file

Command
 ```bash 
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir reference_genome_hg38_index/ --genomeFastaFiles reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --sjdbGTFfile gtf_files/gencode_v39_annotation.gtf --sjdbOverhang 143
 ```
 
Output 
>genome Dir 

 ## Mapping the fastq files to the trimmed reference genome
Input 
>genomeDir, \
>SRR5665260_1.fastp.fastq, \
>SRR5665260_2.fastp.fastq 

Command
```bash
nohup STAR --runThreadN 6 --genomeDir reference_genome_hg38_index/ --readFilesIn SRR5665260_1.fastp.fastq SRR5665260_2.fastp.fastq --sjdbOverhang 143 --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix NA12878_RNA & 
 ```
Output 
>NA12878_RNAAligned.sortedByCoord.out.bam #bam file sorted by coordinate

 ## Mark duplicates
Input 
>NA12878_RNAAligned.sortedByCoord.out.bam 

Command
 ```bash 
 gatk MarkDuplicates --INPUT NA12878_RNAAligned.sortedByCoord.out.bam --OUTPUT NA12878_RNAAligned.sortedByCoord.out_deduped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE NA12878_RNA_duplicates.metrics &
 ```
Output 
>NA12878_RNAAligned.sortedByCoord.out_deduped.bam, \ #deduplicated bam file
>NA12878_RNAAligned.sortedByCoord.out_deduped.bai, \ #deduplicated bam file index
>NA12878_RNA_duplicates.metrics #metrics

 ## SplitNCigar
Input 
>reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta, #reference genome\
>NA12878_RNAAligned.sortedByCoord.out_deduped.bam #deduplicated bam file

Command
 ```bash 
 nohup gatk SplitNCigarReads -R reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I NA12878_RNAAligned.sortedByCoord.out_deduped.bam -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bam &
 ```
Output 
>NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bam, \
>NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bai

 ## Base Recalibration
Input 
>reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta,\
>/workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam,\
>/workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -known-sites,\
>/workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz

Command
 ```bash
 gatk BaseRecalibrator -R /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam --use-original-qualities -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam -known-sites /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -known-sites /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz &
 ```
 
Output 
>NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam 

 ## Apply BQSR
Input 
> /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta,\
> /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam,\
> NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam

Command
```bash
nohup gatk ApplyBQSR --add-output-sam-program-record -R /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam --use-original-qualities -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bam --bqsr-recal-file NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam &
```

Output 
> NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bam,\
> NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bai
 
## Variant calling
 
 



