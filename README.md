# Variant Calling in NA12878 RNA data on AWS
Variant Calling in NA12878 was performed using the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-).

## SETUP A LINUX MACHINE
 ## Sign in to AWS Console (direct sign-in)

### Launch the EC2 dashboard
Search for EC2 in the AWS management console bar
Choose EC2 to open the EC2 Dashboard.

### Create a new key pair
In the left navigation pane, click on Key Pairs. This will display a page to manage your SSH key pairs.

On the Key Pairs page click the Create Key Pair button at the top of the browser window.

In the resulting pop up window, provide a key pair name of your choosing and select .pem
Click on Create key pair.
The key pair you created should automatically download to your system. Follow any browser instructions to save the file to the default download location.
You will see a message appear at the top of the screen that says Successfully created key pair. You will see the key pair you created listed.

##Launch an EC2 instance
Return to the AWS Management Console and open the Amazon EC2 Dashboard.
Click Launch instance, then click Launch instance again from the drop down menu.
Next in the Step 1 page, select the Ubuntu Server 20.04 LTS (HVM), SSD Volume Type - ami-04505e74c0741db8d (64-bit x86) / ami-0b49a4a6e8e22fa16 (64-bit Arm) and click on the Select.
In the Step 2 page, choose an Instance type, select the c5.2xlarge instance from the list and click Next: Configure Instance Details.

On Step 3 page, Configure Instance Details - leave the default settings. Note that the Subnet field can be configured to launch the instance in a specific Availability Zone; while we are keeping the default for this workshop, this gives you control over the location of your machine. Click the Next: Add Storage button in the bottom right corner.

On the Step 4 page, you have the ability to modify or add storage and disk drives to the instance. For this lab, we will simply accept the storage defaults. Take note that the default setting for Delete on Termination is affirmative. This indicates that if the machine is terminated, the root volume associated with the instance will be deleted. You need to uncheck this if you plan to store data on the root volume which you would want to access later. 
When ready, click Next: Configure Security Group.

On Step 6 page, you will be prompted to create a new security group, which will be your firewall rules. Provide a name for your new security group.
Confirm an existing SSH rule exists which allows TCP port 22. To accept connections from Anywhere select the drop-down box under the Source column and select Anywhere which will correspond to 0.0.0.0/0, ::/0.

Click the Review and Launch button.

Review your configuration and choices, and then click Launch.

Select the key pair that you created in the beginning of this lab from the drop-down and check the I acknowledge checkbox. Then click the Launch Instances button.

Your instance will now start, which may take a moment. You will be shown the Launch Status page with the message that your instances are now launching

On the lower right of the page click on View Instances to view the list of EC2 instances. Click on your instance. It will go through an initialization process. Once your instance has launched, you will see your Linux server as well as the Availability Zone the instance is in, and the publicly routable DNS name.

##SSH into EC2 instance
Connecting using SSH on Linux & MacOS and OpenSSH on Windows
In a terminal window, use the ssh command to connect to the instance. Specify the path and file name of the private key (.pem), the user name for your instance, and the public DNS name or IP Address of your instance.

ssh -i /path/my-key-pair.pem ec2-user@<ip-address>
 
 You see a response like the following
 
 The authenticity of host 'ec2-198-51-100-1.compute-1.amazonaws.com (198-51-100-1)' can't be established.
ECDSA key fingerprint is l4UB/neBad9tvkgJf1QZWxheQmR59WgrgzEimCG6kZY.
Are you sure you want to continue connecting (yes/no)?
 
 Enter Yes
You will now be logged into the Instance.
 
 ## Installing conda on AWS
 Using conda to install bioinformatics tools takes care of all dependencies. Users can create different enviroments and install different versions of tools in different environments.
 
 ## Setting up a conda environment on AWS
 Inorder to set up conda on AWS, use the following code.

```bash 
 sudo bash #Give admin rights
 cd / #Go to the root directory
conda create --name=rnaAlign python=3 #This creates a new conda environment with python 3
conda activate rnaAlign #Activate the environment
conda config --add channels bioconda #Add channels
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
In case you would like to expand the partition on the disk use the follwing code.
```bash
growpart /dev/nvme0n1 1
resize2fs /dev/nvme0n1p1
```

 ## Transferring files from local computer to AWS
 

 
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
 
 
 ```bash
scp -i /Users/farheen/Downloads/RNAalign_farheen.pem /Users/farheen/RNA_Alignment_Jan22/genome/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta ubuntu@ec2-3-86-97-53.compute-1.amazonaws.com:/home/ubuntu/
 ```
 
 ## Transferring files from S3 bucket to AWS
 
 
 
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
 
 

 ## Stopping the instance 
You can either go to the Amazon AWS website and search for instance. Select the instance that is running and select Stop instance. This will disconnect you from the instance. 
 
 
 ## Creating a AMI of the instance
 ## Terminating the instance 
 ## Using the AMI to generate an instance 

