#!/bin/bash

for Patient in $(ls /mnt/sda2/TMB/Data/92_patients_NSLC_filtered_VCFS/ | grep ".vcf" | sed 's/NSLC_//' | sed 's/.vcf//'");
do 

    echo "Analysing Patient ${Patient} with VEP"

    vep --af --af_gnomadg --appris --biotype --buffer_size 500 --canonical --ccds --check_existing --distance 5000 --mane --merged --numbers --plugin CADD,/mnt/sda2/VEP_Plugins/CADD/GRCh37/whole_genome_SNVs.tsv.gz,/mnt/sda2/VEP_Plugins/CADD/GRCh37/InDels.tsv.gz --plugin dbscSNV,/mnt/sda2/VEP_Plugins/dbscSNV1.1/dbscSNV1.1_GRCh37.txt.gz --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file /mnt/sda2/TMB/Data/92_patients_NSLC_filtered_VCFS/NSLC_${Patient}.vcf --output_file /mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_${Patient}.vcf --verbose --vcf --port 3337

done
