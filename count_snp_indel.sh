vcftools --gzvcf All_PT_biallelic.vcf.gz --remove indv.txt --recode --stdout | gzip -c > final_dataset_ind.vcf.gz
vcftools --gzvcf All_PT_biallelic.vcf.gz --remove indv.txt --hwe 0.05 --recode --stdout | gzip -c > final_dataset_hwe.vcf.gz
vcftools --gzvcf All_PT_biallelic.vcf.gz --remove indv.txt --minQ 20 --hwe 0.05 --recode --stdout | gzip -c > final_dataset_q20.vcf.gz

gunzip final_dataset_q20.vcf.gz 

vcf2bed --snvs < final_dataset_q20.vcf | wc -l
vcf2bed --insertions < final_dataset_q20.vcf | wc -l 
vcf2bed --deletions < final_dataset_q20.vcf | wc -l
