#!/bin/bash

#USAGE: ./1_process_fastq.pl <plate identifier> 
#Sept 2013 Gregory L Owens

#This pipeline uses BWA/stampy to align and GATK to call SNPs.
#This also demultiplexes 
ID="$1" #A unique identifier for the group of fastq
n_processors="8" #Number of processors to use.

#Files
barcodes="/home/user/SB/Barcodes.$ID.txt" #For GBS demultiplexing. Takes a tab separated file listing sample_name\tbarcode_sequence
ref='/home/user/ref/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa'
stampyref='/home/user/ref/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel'

#Directories
cleandata="/home/user/SB/$ID.clean_data_paired"
trimmeddata="/home/user/SB/$ID.trimmed_data_paired"
unpaired="/home/user/SB/$ID.trimmed_data_unpaired"
sam="/home/user/SB/$ID.sam"
bam="/home/user/SB/$ID.bam"
log="/home/user/SB/$ID.log"
gvcf="/home/user/SB/gvcf"
bin='/home/user/bin'
Hpath='/home/user/SB'
rawdata='/home/user/raw_data'

#Programs
demultiplex='/home/user/bin/1b_GBS_fastq_demultiplexer.pl'
bwa='/home/user/bin/bwa-0.7.9a'
tabix='/home/user/bin/tabix-0.2.6'
vcftools='/home/user/bin/vcftools_0.1.12a/bin'
picardtools='/home/user/bin/picard-tools-1.114'
stampy='/home/user/bin/stampy-1.0.23/stampy.py'
trim="/home/user/bin/Trimmomatic-0.32"
project='SB'
#Prerequisites:
#Index the reference for GATK and BWA
#Change the variables to fit your file structure and ensure the folders exist.

echo "Starting GBS processing pipeline on $ID."
#Make directories if they don't exist
if [ ! -d "$sam" ]; then
	mkdir $sam
fi
if [ ! -d "$bam" ]; then
	mkdir $bam
fi
if [ ! -d "$log" ]; then
	mkdir $log
fi
if [ ! -d "$gvcf" ]; then
        mkdir $gvcf
fi
if [ ! -d "$cleandata" ]; then
        mkdir $cleandata
fi
if [ ! -d "$trimmeddata" ]; then
        mkdir $trimmeddata
fi
if [ ! -d "$unpaired" ]; then
        mkdir $unpaired
fi

#Demultiplex GBS data. Comment out if data is already demultiplexed and in $cleandata folder.
perl $demultiplex $barcodes $rawdata/"$ID"_R1.fastq $rawdata/"$ID"_R2.fastq $cleandata/

#Create list of samples.
ls $cleandata | grep -v "nobar" | sed s/_R1//g | sed s/_R2//g | sed s/.fastq// | uniq  > $Hpath/Samplelist.$ID.txt

#Trim the data using Trimmomatic. Removes bad reads and illumina adapter contamination.
while read prefix
do
java -jar $trim/trimmomatic-0.32.jar PE -phred33 $cleandata/"$prefix"_R1.fastq $cleandata/"$prefix"_R2.fastq $trimmeddata/"$prefix"_R1.fastq $unpaired/"$prefix"_unR1.fastq $trimmeddata/"$prefix"_R2.fastq $unpaired/"$prefix"_unR2.fastq ILLUMINACLIP:$trim/adapters/TruSeq3-PE.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36
done < Samplelist.$ID.txt


###Align using BWA. Turn from sam to bam. Sort by coordinate and add read group data.
while read prefix
do
	#Align paired reads using bwa aln and sampe
	$bwa/bwa aln -t $n_processors $ref $trimmeddata/"$prefix"_R1.fastq 1> $trimmeddata/"$prefix"_R1.sai
	$bwa/bwa aln -t $n_processors $ref $trimmeddata/"$prefix"_R2.fastq 1> $trimmeddata/"$prefix"_R2.sai
	$bwa/bwa sampe $ref $trimmeddata/"$prefix"_R1.sai $trimmeddata/"$prefix"_R2.sai $trimmeddata/"$prefix"_R1.fastq $trimmeddata/"$prefix"_R2.fastq 1> $sam/$prefix.sam 2> $log/$prefix.bwasampe.log
	#Convert to bam
	samtools view -Sb $sam/$prefix.sam > $bam/${prefix}_premerge.bam
	#Realign using stampy
        $stampy -g $stampyref -h $stampyref -t$n_processors --bamkeepgoodreads -M $bam/${prefix}_premerge.bam -o $bam/${prefix}_premerge.stampy.bam 2> $log/$prefix.stampy.log
	
	#Align unpaired read 1 using bwa aln and samse
	$bwa/bwa aln -t $n_processors $ref $unpaired/"$prefix"_unR1.fastq 1> $unpaired/"$prefix"_unR1.sai
	$bwa/bwa samse $ref $unpaired/"$prefix"_unR1.sai $unpaired/"$prefix"_unR1.fastq 1> $sam/${prefix}_unR1.sam 2> $log/${prefix}_unR1.bwasamse.log
	#Convert to bam
	samtools view -Sb $sam/${prefix}_unR1.sam > $bam/${prefix}_unR1.bam
	#Realign using stampy
	$stampy -g $stampyref -h $stampyref -t$n_processors --bamkeepgoodreads -M $bam/${prefix}_unR1.bam -o $bam/${prefix}_unR1.stampy.bam 2> $log/${prefix}_unR1.stampy.log	
	
	#Align unpaired read 2 using bwa aln and samse
	$bwa/bwa aln -t $n_processors $ref $unpaired/"$prefix"_unR2.fastq 1> $unpaired/"$prefix"_unR2.sai
	$bwa/bwa samse $ref $unpaired/"$prefix"_unR2.sai $unpaired/"$prefix"_unR2.fastq 1> $sam/${prefix}_unR2.sam 2> $log/${prefix}_unR2.bwasamse.log
	#Convert to bam
	samtools view -Sb $sam/${prefix}_unR2.sam > $bam/${prefix}_unR2.bam
	#Realign using stampy
	$stampy -g $stampyref -h $stampyref -t$n_processors --bamkeepgoodreads -M $bam/${prefix}_unR2.bam -o $bam/${prefix}_unR2.stampy.bam 2> $log/${prefix}_unR2.stampy.log
	
	#Merge bam files for all 3 sets
	java -jar $picardtools/MergeSamFiles.jar INPUT=$bam/${prefix}_premerge.stampy.bam INPUT=$bam/${prefix}_unR1.stampy.bam INPUT=$bam/${prefix}_unR2.stampy.bam OUTPUT=$bam/${prefix}.sort.bam SORT_ORDER=coordinate USE_THREADING=true 2> $log/$prefix.mergesam.log 
	#Clean bam files
	java -jar $picardtools/CleanSam.jar INPUT=$bam/$prefix.sort.bam OUTPUT=$bam/$prefix.clean.bam 2> $log/$prefix.cleansam.log
	#Add read group information (needed for GATK)
	java -jar $picardtools/AddOrReplaceReadGroups.jar I=$bam/$prefix.clean.bam O= $bam/$prefix.sortrg.bam SORT_ORDER=coordinate RGID=$prefix RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$prefix CREATE_INDEX=True 2> $log/$prefix.addRG.log

	#Remove intermediate files
	rm $trimmeddata/"$prefix"_R1.sai
	rm $trimmeddata/"$prefix"_R2.sai
	rm $unpaired/"$prefix"_unR1.sai
	rm $unpaired/"$prefix"_unR2.sai
	rm $sam/$prefix.sam
	rm $sam/${prefix}_unR1.sam
	rm $sam/${prefix}_unR2.sam
	rm $bam/${prefix}_premerge.bam
	rm $bam/${prefix}_unR1.bam
	rm $bam/${prefix}_unR2.bam
	rm $bam/${prefix}_premerge.stampy.bam
	rm $bam/$prefix.stampy.bam
	rm $bam/${prefix}_unR1.stampy.bam
	rm $bam/${prefix}_unR2.stampy.bam
	rm $bam/$prefix.clean.bam
	rm $bam/$prefix.sort.bam
done < Samplelist.$ID.txt


#Make bam.list for GATK
ls -d $bam/*.* | grep sortrg.bam  > $Hpath/bamlist.$ID.list

#identify local indels

java -Xmx4g -jar $bin/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R $ref \
   -I $Hpath/bamlist.$ID.list \
   -nt $n_processors \
   -log $log/$ID.RealignerTargetCreator.log \
   -o $Hpath/$ID.realign.intervals

#Realign around local indels
while read prefix
do
java -Xmx4g -jar $bin/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $ref \
	-I $bam/$prefix.sortrg.bam \
	-targetIntervals $Hpath/$ID.realign.intervals \
	-o $bam/$prefix.realign.bam \
	-log $log/$ID.$prefix.IndelRealigner.log 

#Call GATK HaplotypeCaller to make gvcfs
java -Xmx18g -jar $bin/GenomeAnalysisTK.jar \
	-nct $n_processors \
	-l INFO \
	-R $ref \
	-log $log/$ID.$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-I  $bam/$prefix.realign.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 2 \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o $gvcf/${prefix}.GATK.gvcf.vcf
done < $Hpath/Samplelist.${plate}.txt

#Make input list for GATK GenotypeGVCFs
tmp=""
while read prefix
do
        tmp="$tmp --variant $gvcf/$prefix.GATK.gvcf.vcf"
done < $Hpath/Samplelist.$ID.txt

#Genotype all gvcf together into one vcf file
java -Xmx18g -jar $bin/GenomeAnalysisTK.jar \
	-nt $n_processors \
	-l INFO \
	-R $ref \
	-log $log/$ID.GenotypeGVCFs.log \
	-T GenotypeGVCFs \
	$tmp \
	-o $Hpath/SB.GATK.total.vcf \
	-inv \
	--max_alternate_alleles 4


exit

