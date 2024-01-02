# Ancestral State Inference
Adapted from script

## 1. Process input files, the sampling populations and the ancestral group population. 
To cut down the running time, I splitted the files into 58 chromosomal regions with job-array. See script.sh for detail.
If two groups were called together, this step can be skipped. 

1.a) To remove AD information for both ```sample_group.vcf``` and ```anc_group.vcf```. These VCFs need to include both variants and invariants (quality filtered
```
bcftools annotate -x FORMAT/AD sample_group.vcf \
          -O z \
          -o sample_group.noAD.vcf.gz

tabix sample_group.noAD.vcf.gz
```

1.b) To merge the VCFs

```
bcftools merge sample_group.noAD.vcf.gz anc_group.noAD.vcf.gz \
          -O z \
          -o 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz

tabix 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz
```
#Remove multiallelic sites

```
vcftools --gzvcf 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz \
        --min-alleles 2 --max-alleles 2 \
        --max-missing-count 105 \
        --recode --out 116combineVCFs/116.${chrm}_${startpos}_${endpos}.biallelic.missing_removed
bgzip 116combineVCFs/116.${chrm}_${startpos}_${endpos}.biallelic.missing_removed.recode.vcf
tabix 116combineVCFs/116.${chrm}_${startpos}_${endpos}.biallelic.missing_removed.recode.vcf.gz

```
## 2. Use python script to infer the ancestral state. 
```
python3 anc_state_infer.py \
        -I 116combineVCFs/116.${chrm}_${startpos}_${endpos}.biallelic.missing_removed.recode.vcf.gz \
        -O 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf \
        -Alist GJF.popline.txt
```
