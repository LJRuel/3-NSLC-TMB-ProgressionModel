#!/bin/bash

#@author Louis-Jacques Ruel, louis.jacques.ruel@gmail.com 01/2023

#This script runs VEP on a series of already annotated vcf files. This script add the --pick option of VEP if not included in the first VEP annotation. This results in only one annotation per variant (instead of multiple for each variant).
#Requires VEP to be installed with CADD and dbscSNV plugins.


# Change vcf list (here all patients ID)
for Patient in $(cat /mnt/sda2/TMB/Data/patients_ID_list.txt);
do 

    printf "\nAdjusting VEP annotated variants Patient ${Patient} with VEP.\n"

    # Change the --pick order if needed, the input and output. See https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html and http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick for --pick documentation.

    vep -i /mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_${Patient}.vcf -o /mnt/sda2/TMB/Data/VEP_pick_92_NSLC/NSLC_${Patient}.vcf --cache --merged --pick --force_overwrite --vcf --species homo_sapiens --verbose --port 3337

done

