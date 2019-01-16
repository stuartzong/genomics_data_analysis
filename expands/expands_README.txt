1. Introduction:
The Expands pipeline does the following:
-> Combine DNA and RNA variants for tumor-normal pairs from mpileup, mutseq, strelka pipeline etc. 
We look for variants in coding regions and ignore ones in intron and intergenic regions. Theoreticall, it also makes sense to include variants classfied as MODIFIER by SNPEFF. Although it is unlikely modifier variants will provide selective advantage, they could still be used as markers to distinguish different cell populations. 

-> Summarize variants by gene and variant level to create a temporary summary file. 
The python script parses all vcf files including DNA and RNA relevant to a patient and summarize the variants in both gene and variant level. This is useful when we are interested in looking at the recurrence of mutations. 

-> Rerun samtools mpileup at default setting for all variant positions.
For each patient and tissue status, the python script generates a mpileup bash script, which is automatically submitted to m0001 cluster. The script repeatly detect if the cluster jobs are finished by checking the completion stamps.
  
-> Parse mpileup results to get accurate base count info for reference and alternative base, calculte allele frequency.
The python script invokes another python script to parse the mpileup raw output files and report the coverage, reference base count, and alternative base count for each variant position.
 
-> Combine base count info for tumour-normal pairs if applicable.
If there is matched normal for a tumor, the script provide the the coverage, reference base count, and alternative base count for both tumor and normal.

-> Adjust tumour base count according to tumour content if provided.
If tumor content is known, the script adjusts the allele frequecy  based on the tumor content. 
 
-> Perform Fisher Exact test to determine p value.
The script calls a R script to calculate the Fisher Exact Test p value which indicates the allele frequency difference between tumor and normal is statistically significant.

-> Annotate the temprary summary file to create a final SNV summary file.
This incorporate the base counts and the caller(s) which has called the variant into the temporary summary file.
 
-> Perfom confidence and quality based filtering to generate a list of high confident somatic mutations. 
It is important to only include high confidence variant (for example af > 10%) because expands is intended for low to moderate coverage genome data. Within this context, variants eith lower than 10% allele frequence are very often false positive. Also, to achieve stable result, a minimum of 200 variants is recommended. 

-> Generate expands input files.
Expands takes in two input files: somatic SNVs and copy numbers.  SNV file is made by reformatting  the SNV summary and copy number file is made by copying the CNAseq seg.txt file. Please note expands takes subclonal copy number. Therefore, more accurate copy number estimate is prefered. But here we use CNAseq integer estimates instead since this is the only data available.

-> Generate R scripts to run expands, which estimates cellular frequency and cluster variants into subpopulations.

2. EXPANDS: Expanding Ploidy and Allele-Frequency on Nested Subpopulations
Expanding Ploidy and Allele Frequency on Nested Subpopulations (expands) characterizes coexisting subpopulations in a single tumor sample using copy number and allele frequencies derived from exome- or whole genome sequencing input data (http://www.ncbi.nlm.nih.gov/pubmed/24177718). The model detects coexisting genotypes by leveraging run-specific tradeoffs between depth of coverage and breadth of coverage. This package predicts the number of clonal expansions, the size of the resulting subpopulations in the tumor bulk, the mutations specific to each subpopulation and tumor purity. The main function runExPANdS() provides the complete functionality needed to predict coexisting subpopulations from single nucleotide variations (SNVs) and associated copy numbers. The robustness of subpopulation predictions increases with the number of mutations provided. It is recommended that at least 200 mutations are used as input to obtain stable results. Updates in version 1.6 include: (1) So far mutations had been assigned to maximal one subpopulation. However mutations may not be exclusive to the assigned subpopulation but may also be present in smaller, descending subpopulations. Whether or not this is the case is now decided by leveraging the predicted phylogenetic structure of the subpopulation composition. (2) Included homozygous deletion as potential scenario when modeling (SNV,CNV) pairs with overlapping genomic location, that are propagated during distinct clonal expansions. (3) Optimized solution to improve sensitivity at cell-frequency distribution margins. Need for improvement was because subpopulation detection sensitivity correlates to centrality of subpopulation size during clustering. Tolerance of copy number and allele frequency measurement errors must be higher for marginal cell-frequencies than for central cell-frequencies, in order to counteract the reduced cluster detection sensitivity at the cell-frequency distribution margins. This is only relevant during subpopulation detection (SNV clustering), cell-frequency independent error tolerance still applies during SNV assignment. (4) Fixed a bug where incorrect data matrix conversion could occur when handing non-numerical matrix as parameter to function runExPANdS(). Further documentation and FAQ around this package is available at http://dna-discovery.stanford.edu/software/expands.

3. Pipeline input files: bam_vcf_cnv_path.txt (tab-delimited).
Information such as patient name, tissue status, and full file paths to bam, copy number segmentation, single or paired mpielup vcf, mutationseq vcf, and strelka vcf, and tumor content etc must be specified in the input file. The input file must contain these headers: patient, DNA_lib, RNA_lib, status, tumour_content, bam_path, paired_mpileup_vcf, single_vcf, mutseq_snv_vcf, strelka_snv_vcf, strelka_indel_vcf, cnv, RNA_bam, RNA_single_vcf, other_vcf.

The python scripts generate expands input files: SNV matrix and copy number matrix:
-> SNV matrix: 
It is recommended to include minimum 200 SNVs to achieve stable results. By default, 8000 SNVs are the maximum. If more than 8000 SNVs are provided, expands will only cluster SNVs in the SureSelectExome_hg19, comprising ca. 468 MB centered on the human exome. Also, only autosomal SNVs will be accounted for. Each row corresponds to a point mutation. Only mutations located on autosomes should be included. Columns in SNV must be labeled and must include:
chr- the chromosome on which each mutation is located;
startpos- the genomic position of each mutation;
AF_Tumor- the allele-frequency of each mutation;
PN_B- ploidy of B-allele in normal cells. A value of 0 indicates that the mutation has only been detected in the tumor sample (i.e. somatic mutation). A value of 1 indicates that the mutation is also present in the normal (control) sample, albeit at reduced allele frequency (i.e. mutation is consequence of LOH). Mutations, for which the allele frequency in the tumor sample is lower than the corresponding allele frequency in the normal sample, should not be included.

-> copy number matrix:
Each row corresponds to a copy number segment. Copy number matrix is typically the output of a circular binary segmentation algorithm. Columns in CBS must be labeled and must include:
chr- chromosome;
startpos- the first genomic position of a copy number segment;
endpos- the last genomic position of a copy number segment;
CN_Estimate- the absolute copy number estimated for each segment

4. Run expands R script
Install expands R package and run the R scripts to generate expands graphs and subpopulations.

5. Result interpretation
For each point mutation (x-axis) the function displays:
-> the size of the subpopulation to which the mutation has been assigned (squares). Each square is colored based on the confidence with which the mutation has been assigned to the corresponding subpopulation (black - highest, white - lowest).
-> the total ploidy of all alleles at the mutated genomic locus in that subpopulation (dots).
-> the allele frequency of the mutation. Allele frequencies and ploidities are  colored based on the chromosome on which the mutation is located (stars - somatic mutations, triangles - loss of heterozygosity).
