####! /bin/bash
####$ -S /bin/bash
####$ -N vcf2maf_{{patient}}
####$ -q transabyss.q
####$ -l mem_token=2G,mem_free=2G,h_vmem=2G
####$ -V
###
#source /gsc/software/linux-x86_64-centos6/vcf2maf-1.6.6/profile.sh
mkdir /projects/trans_scratch/validations/workspace/szong/Cervical/mutsigcv/mafs/{{patient}}
cd /projects/trans_scratch/validations/workspace/szong/Cervical/mutsigcv/mafs/{{patient}} 

/gsc/software/linux-x86_64-centos6/perl-5.22.0/bin/perl /gsc/software/linux-x86_64-centos6/vcf2maf-1.6.6/vcf2maf.pl --input-vcf  {{input_vcf}} --output-maf {{output_maf}} --tumor-id  {{patient}}_T --normal-id {{patient}}_N --vep-path /gsc/software/linux-x86_64-centos6/vcf2maf-1.6.6/vep --vep-data /projects/analysis/vepdata --vcf-tumor-id TUMOR --vcf-normal-id NORMAL --custom-enst /gsc/software/linux-x86_64-centos6/vcf2maf-1.6.6/data/isoform_overrides_uniprot --output-directory /projects/trans_scratch/validations/workspace/szong/Cervical/mutsigcv/mafs/{{patient}}

if [ $? -eq 0 ]
then
    echo "vcf2maf job finished successfully at:" `date`
else
    echo "ERROR: vcf2maf did not finish correctly!"
    exit 3
fi

