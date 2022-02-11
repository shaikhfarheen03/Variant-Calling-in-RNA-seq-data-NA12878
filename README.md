# Variant Calling in NA12878 RNA data on AWS
Variant Calling in NA12878 was performed using the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-).

# Fastq files and reference genome
Fastq files were obtained from the following link https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5665260 
To download the files I installed SRAtoolkit on the server and used the commmand fasterq dump command: 
```bash
fasterq-dump SRR5665260
```

fasterq-dump also takes --threads and --progress as parameters. 

Reference genome was downloaded on the local computer from the following link: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false 
The fasta file, .fai file, .dict, dbsnp file, known indel sites were downloaded.

GTF file was downloaded from Gencode: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz

 
 # SETUP A LINUX MACHINE
 ## Sign in to AWS Console (direct sign-in)
![Screen Shot 2022-02-11 at 3 29 07 PM](https://user-images.githubusercontent.com/26681884/153665861-29e792e4-403d-438f-94db-655e9132a0f0.png)

## Launch the EC2 instance
Search for EC2 in the AWS management console bar
![Screen Shot 2022-02-11 at 3 33 58 PM](https://user-images.githubusercontent.com/26681884/153666386-a27e5984-1c3a-4963-a535-a040f803e062.png)

Choose EC2 to open the EC2 Dashboard.

 
 
 
 
 After creating an account generate a .PEM key by going to 
 The key will be downloaded on your local computer
 
 ## Starting a AWS EC2 instance
 Start an EC2 instance and add enough memory to the instance. 
 
 ## Logging in to the instance 
 Go to terminal and type 
 
 
 ## Installing conda on AWS
 ## Setting up a conda environment on AWS
```bash 
conda create --name=rnaAlign python=3
conda activate rnaAlign
conda config --add channels bioconda
conda config --add channels conda-forge
```

 ## Installing common bioinformtic tools
```bash
conda install samtools
conda install gatk
conda install fastp
conda install STAR
conda install 
```

 ## Expanding partition on disk
 
```bash
growpart /dev/nvme0n1 1
resize2fs /dev/nvme0n1p1
```

 ## Transferring files from local computer to AWS
 
```bash
scp -i /Users/farheen/Downloads/RNAalign_farheen.pem /Users/farheen/RNA_Alignment_Jan22/genome/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta ubuntu@ec2-3-86-97-53.compute-1.amazonaws.com:/home/ubuntu/
 ```
 
 ## Transferring files from S3 bucket to AWS
 
 
 
 ## Stopping the instance 
You can either go to the Amazon AWS website and search for instance. Select the instance that is running and select Stop instance. This will disconnect you from the instance. 
 
 
 ## Creating a AMI of the instance
 ## Terminating the instance 
 ## Using the AMI to generate an instance 
 # Workflow
 ### Running fastp for qaulity and to trim adapters (You could also use fastqc to get quality and then clean the fastq using fastp)
 - Input SRR5665260_1.fastq, SRR5665260_2.fastq
```bash 
fastp -i SRR5665260_1.fastq -I SRR5665260_2.fastq -o SRR5665260_1.fastp.fastq -O SRR5665260_2.fastp.fastq --detect_adapter_for_pe 
```
- Output SRR5665260_1.fastp.fastq, SRR5665260_2.fastp.fastq, fastp.html, fastp.json
- Insert fastp report

 ### Generating a index for the reference genome
 Why creating a index for the reference genome is important
  - Input resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta, gencode_v39_annotation.gtf
 ```bash 
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir reference_genome_hg38_index/ --genomeFastaFiles reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --sjdbGTFfile gtf_files/gencode_v39_annotation.gtf --sjdbOverhang 143
 ```
- Output genome Dir 

 ### Mapping the fastq files to the trimmed reference genome
 Some theory about STAR and 2 pass basic mode

- genome index - Input genomeDir, SRR5665260_1.fastp.fastq, SRR5665260_2.fastp.fastq 
```bash
nohup STAR --runThreadN 6 --genomeDir reference_genome_hg38_index/ --readFilesIn SRR5665260_1.fastp.fastq SRR5665260_2.fastp.fastq --sjdbOverhang 143 --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix NA12878_RNA & 
 ```
- Output NA12878_RNAAligned.sortedByCoord.out.bam

 ### Mark duplicates
 Some theory about Mark Duplicates
 - Input NA12878_RNAAligned.sortedByCoord.out.bam
 ```bash 
 gatk MarkDuplicates --INPUT NA12878_RNAAligned.sortedByCoord.out.bam --OUTPUT NA12878_RNAAligned.sortedByCoord.out_deduped.bam --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --METRICS_FILE NA12878_RNA_duplicates.metrics &
 ```
- Output NA12878_RNAAligned.sortedByCoord.out_deduped.bam, NA12878_RNAAligned.sortedByCoord.out_deduped.bai, NA12878_RNA_duplicates.metrics

 ### SplitNCigar
 Some theory about Split and Cigar
  - Input reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta, NA12878_RNAAligned.sortedByCoord.out_deduped.bam
 ```bash 
 nohup gatk SplitNCigarReads -R reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I NA12878_RNAAligned.sortedByCoord.out_deduped.bam -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bam &
 ```
- Output NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bam, NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar.bai

 ### Base Recalibration
 Some theory about Base Recalibration
 - Input reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta, /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam, /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -known-sites, /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz
 ```bash
 gatk BaseRecalibrator -R /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam --use-original-qualities -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam -known-sites /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz -known-sites /workspace/RNAseq_data/gatk_resource_bundle/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz &
 ```
- Output NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam 

 ### Apply BQSR
 Some theory about BQSR 
 - Input /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta, /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam, NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam
```bash
nohup gatk ApplyBQSR --add-output-sam-program-record -R /workspace/RNAseq_data/reference_genome_hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta -I /workspace/RNAseq_data/NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp.bam --use-original-qualities -O NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bam --bqsr-recal-file NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal.bam &
```
- Output NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bam, NA12878_RNAAligned.sortedByCoord.out_deduped_postSplitNCigar_withReadGrp_recal_2.bai
 
# Variant calling
 
 


