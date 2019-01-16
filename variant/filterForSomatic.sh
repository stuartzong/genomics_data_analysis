#!/bin/bash -x

ID=$(basename $1)

#Get just the coding variants from the initial file
grep MISSENSE $1 >  $ID.missense.vcf
grep NONSENSE $1 >  $ID.nonsense.vcf
grep ^# $1 > $ID.header
cat $ID.missense.vcf $ID.nonsense.vcf | sort -u | cat $ID.header - |/home/rcorbett/bin/vcftools_0.1.11/bin/vcf-sort -c > $ID.coding.vcf

/home/rcorbett/slx_service_rc/trunk/germlineSomatic.py -v $ID.coding.vcf -r /projects/rcorbettprj2/germline_or_somatic/analysis/vcfs_fields.txt > $ID.coding.txt

#Keep just the SNVs (INDELS will have to be handled differently - sigh) - reformat for matching
awk '$1~/[0-9A-Z]+_[0-9]+_[ACTG]_[ACTG]$/ { printf "%s %s %s %s %s %s %s\n", $1,$2,$3,$4,$5,$6,$7 }' $ID.coding.txt > $ID.coding.justSNPs.txt

#Our best efficacy comes from these two lines
awk '$2!="-" || $3!="-" || $4!="-"' $ID.coding.justSNPs.txt > $ID.coding.justSNPs.notedGermline.txt
awk '$2=="-" && $3=="-" && $4=="-"' $ID.coding.justSNPs.txt > $ID.coding.justSNPs.notGermline.txt 
awk '$5!="-" || $6!="-" || $7!="-"' $ID.coding.justSNPs.txt > $ID.coding.justSNPs.notedSomatic.txt #Maybe include common somatic variants

#Keepin the lines below in case we need to go back to these in the future.
#awk '$5!="-" || $6!="-" || $7!="-"' $ID.notGermline.txt > $ID.notGermline.notedSomatic.txt
#awk '$5=="-" && $6=="-" && $7=="-"' $ID.notGermline.txt > $ID.notGermline.notSomatic.txt
#awk '$5!="-" || $6!="-" || $7!="-"' $ID.notedGermline.txt > $ID.notedGermline.notedSomatic.txt
#awk '$5=="-" && $6=="-" && $7=="-"' $ID.notedGermline.txt > $ID.notedGermline.notSomatic.txt
#awk '$5!="-" || $6!="-" || $7!="-"' $ID.justSNPs.txt > $ID.notedSomatic.txt
#awk '$5=="-" && $6=="-" && $7=="-"' $ID.justSNPs.txt > $ID.notSomatic.txt

grep ^# $1 > $ID.header  #Keep the original VCF header so we can put it back at the top of the result

# Collect the germline matched variants in ExAC  with AF<=0.0001 and put them together with the variants that aren't in any of the germline references.
# There's got to be a cleaner way to do this
awk '{ print $1, $4, $8 }' $ID.coding.justSNPs.notedGermline.txt | sed 's/AF=//' | awk '$2!="-"' | awk '$2 < 0.0001' | cat $ID.coding.justSNPs.notGermline.txt - | tr '_' '\t' | awk '{ print $1, $2 }' | tr ' ' '\t' | xargs -i grep -w "{}" $ID.coding.vcf | cat $ID.header - | /home/rcorbett/bin/vcftools_0.1.11/bin/vcf-sort -c > $ID.maybeSomaticSNVS.vcf

