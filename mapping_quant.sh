# Description:
# This script runs the following steps:
# 1. Run fastQC to QC files
# 2. Run STAR to map reads to genome/transcriptome - HUMAN index, genecode GRCh38
# 3. Run featureCounts to generate count table - Human GTF
# 4. run multiQC to summarize results

# specify
group_project=arpe19_miom1 #project name

star_index=/data/reference/human/star_index_100bp #STAR mapping index

GTF=/data/reference/human/gencode.v30.primary_assembly.annotation.gtf

fastqc_threads=20
cutadapt_threads=20
star_threads=5
fc_threads=20

# general setup
fastq_path=/data/primary/${group_project}
output_path=/data/secondary/${group_project}

mkdir -p $output_path

# fastqc, run for each fastq separately
mkdir $output_path/fastqc
fastqc -o $output_path/fastqc -t $fastqc_threads $fastq_path/*.fastq

# multiqc summary of fastqc html
mkdir $output_path/multiqc
mkdir $output_path/multiqc/fastqc
multiqc -f -o $output_path/multiqc/fastqc $output_path/fastqc


# STAR with human genome index
echo "Start STAR"
mkdir $output_path/star

for fastq_i in `ls ${fastq_path}/*.fastq`
do
	filename=`basename $fastq_i`
	filename=${filename%.fastq}_
	
	STAR --runThreadN $star_threads \
	--genomeDir $star_index \
	--outFileNamePrefix $output_path/star/$filename \
	--outSAMtype BAM SortedByCoordinate \
	--readFilesIn $fastq_i
done

# multiqc summary
mkdir $output_path/multiqc/star
multiqc -f -o $output_path/multiqc/star $output_path/star


# Run featureCounts to quantify gene expression of STAR aligned bam
echo "Start featureCounts"
mkdir $output_path/star_featureCounts

for bam_i in `ls $output_path/star/*.bam`
do
	filename=`basename $bam_i`
	filename=${filename%_Aligned.sortedByCoord.out.bam}
	featureCounts -T $fc_threads \
	-a $GTF \
	-o $output_path/star_featureCounts/$filename.counts $bam_i
done

# multiqc summary
mkdir $output_path/multiqc/featureCounts
multiqc -f -o $output_path/multiqc/featureCounts $output_path/star_featureCounts

