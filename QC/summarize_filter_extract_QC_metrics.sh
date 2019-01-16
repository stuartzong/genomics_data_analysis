#!usr/bin/sh
#This script retrieves bam to genome bam files for RNAseq libraries
#This script is used to create input file for RNAseq QC analysis from bam file list
#parse all QC files to generate QC report

echo "Please run this script on xhost!"

wkdir=`pwd`

libs=$1
project=$2
bam="bam_files.txt"
flowcell="lib.flowcell.pool.txt"
qcinput="rnaqc.input.txt"
script_path="/home/szong/coverage-based-QC/bin"

echo  "get bam path for a list of libraries."
rm $bam

while read line; do
    bam_path=$(bash /home/rmar/bin/retrieve_reads_to_genome_bam_file $line|grep "withJunctionsOnGenome_dupsFlagged.bam$"|sed 's/hg19 BAM file found: //g'|sed 's/BAM file found: //g')
    if [[ ! -z ${bam_path//[[:blank:]]/} ]]
       then
       echo $bam_path >> $bam
    fi
done < $libs

awk -F "/" '{print $7"\t"$6"\t"$5}' $bam > $flowcell


echo "Making rnaqc input files."
rm $qcinput
while read line; do
  echo $line
  lib=$(echo $line|awk '{print $1}')
  pool=$(echo $line|awk '{print $3}'|tr '[:upper:]' '[:lower:]')
  flowcell=$(echo $line|awk '{print $2}'|tr '[:upper:]' '[:lower:]')
  echo "the library id is: $lib"
  #coverage_QC_path="/projects/wtsspipeline/analysis/$lib/coverage/hg19a_jg-e69/ensembl_homo_sapiens_core_69_37/production_v1.1/no_strand/"
  # bwa_mem used since January 2016
  coverage_QC_path="/projects/wtsspipeline/analysis/$lib/coverage/hg19a_jg-e69_mem/ensembl_homo_sapiens_core_69_37/production_v1.1/no_strand/"
  stranded_coverage_QC_path="/projects/wtsspipeline/analysis/$lib/coverage/hg19a_jg-e69_mem/ensembl_homo_sapiens_core_69_37/production_v1.1/no_strand/"
  #snp_path="/projects/wtsspipeline/analysis/$lib/external/samtools_mpileup/hg19a_jg-e69/mafs_v1.1/rna_tumor/"
  snp_path=""
  p5_p3_ratio_path="/projects/wtsspipeline/analysis/$lib/5-3_ratio_revision99393/hg19a_jg-e69_mem/ensembl_homo_sapiens_core_69_37/production_v1.1/"
  #level0_qc_path="/projects/tab/QCreports/ix2189/c34ucacxx_5/$lib/"
  level0_qc_path="/projects/tab/QCreports/$pool/$flowcell/$lib/"
  leakage_path=""
  (echo $lib","$coverage_QC_path","$stranded_coverage_QC_path","$snp_path","$p5_p3_ratio_path","$level0_qc_path","$leakage_path",") >> $qcinput
done < $flowcell

echo "run python scripts to gather all QC metrics from input files"
echo "Generating RNA QC report!"
cd $script_path
source /home/szong/coverage-based-QC/environment/coverage_based_qc_dev_env.sh
python /home/szong/coverage-based-QC/bin/level1_qc_report_main.py -c $wkdir/$qcinput -o $wkdir -f -p $project
cd $wkdir

# extract fields from qc report
#1	library_name
#5	percent_aligned_of_chastity_passed
#3	total_gbp_chastity_passed
#20	exon_intron_ratio
#17	num_genes_coverage_ge_10X

echo "Parsing QC metrices to extract info including: library_name, Aligned_Gbp, total_gbp, exon_intron_ratio, and num_genes_coverage_ge_10X!"
date=`date +%Y%m%d`
echo $date
report=$project"_level1_qc_summary_report_"$date".verbose.txt"

awk -F "\t" '{print $2"\t"$28"\t"$23"\t"$48"\t"$45}' $report > $report".final"
