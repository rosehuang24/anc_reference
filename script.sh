#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p amd-ep2
#SBATCH -J split107
#SBATCH -o ./report/split107.%A_%a.out
#SBATCH -e ./report/split107.%A_%a.error
#SBATCH -a 1-58

source ~/.bash_profile

module unload gcc
module load gcc/9.3.0
module load jdk/1.8.0
module load fastqc/0.11.9
module load trimmomatic/0.39
module load bwa/0.7.17
module load gatk/4.2.0.0
module load picard/2.25.1
module load perl/5.34.0
module load samtools/1.13
module load sratoolkit/2.11.2
module load htslib/1.12
module load vcftools/0.1.16
module load bcftools/1.14

parentDIR=/storage/zhenyingLab/huangruoshi
CRNTDIR=$parentDIR/anc_inferrence
DATADIR=$parentDIR/CALLs_scattered/GJF/filtered_vcf
TXTDIR=/storage/zhenyingLab/huangruoshi/txt_might_be_useful
REFDIR=/storage/zhenyingLab/huangruoshi/chicken_ref
#INDIR=$CRNTDIR/raw_output_vcf
#OUTDIR=$CRNTDIR/filtered_vcf

ref_genome=$REFDIR/Gallus_gallus.GRCg6a.dna.toplevel.fa
ref_bwa=$REFDIR/RefSeq
dbsnp=$REFDIR/newdbSNP.vcf.gz

regionfile=$TXTDIR/GGchrm-regions.txt
#chrm=`head -n ${SLURM_ARRAY_TASK_ID} $regionfile |tail -n1 | awk '{print $1}'`
#region=`head -n ${SLURM_ARRAY_TASK_ID} $regionfile |tail -n1 | awk '{print $2}'`
regionname=`head -n ${SLURM_ARRAY_TASK_ID} $regionfile |tail -n1 | awk '{print $3}'`

coordfile=$TXTDIR/GG_coords.txt
chrm=`head -n ${SLURM_ARRAY_TASK_ID} $coordfile | tail -n1 | awk '{print $1}'`
startpos=`head -n ${SLURM_ARRAY_TASK_ID} $coordfile | tail -n1 | awk '{print $2}'`
endpos=`head -n ${SLURM_ARRAY_TASK_ID} $coordfile | tail -n1 | awk '{print $3}'`
#



SNPVCF=$parentDIR/108/snp.input.4.masked_removed.recode.vcf.gz

vcftools --gzvcf $SNPVCF \
       --chr $chrm --from-bp $startpos --to-bp $endpos \
       --recode \
       --out tmp/input.4.107.${chrm}_${startpos}_${endpos}

bcftools annotate -x FORMAT/AD tmp/input.4.107.${chrm}_${startpos}_${endpos}.recode.vcf -O z -o tmp/input.4.107.noAD.${chrm}_${startpos}_${endpos}.vcf.gz
poolAD=tmp/input.4.107.noAD.${chrm}_${startpos}_${endpos}.vcf.gz
tabix $poolAD

#those are already with singletons. further pipemark-process is futile. 
gjf=$parentDIR/CALLs_scattered/GJF/filtered_vcf/snp_hom_allsites_${regionname}.vcf.gz

bcftools annotate -x FORMAT/AD $gjf -O z -o tmp/gjfallsites_filtered.noAD.${chrm}_${startpos}_${endpos}.vcf.gz
gjfAD=tmp/gjfallsites_filtered.noAD.${chrm}_${startpos}_${endpos}.vcf.gz
tabix $gjfAD

#you can actually do it in earlier steps.
vcftools --gzvcf $gjfAD \
        --exclude-bed $REFDIR/galGal6.masked.bed \
        --recode \
        --out tmp/gjfallsites_repRemoved.noAD.${chrm}_${startpos}_${endpos}
        
##NOTE!!!: this maskRemoved vcf still have singleton. Take care of it in python script!       
bgzip tmp/gjfallsites_repRemoved.noAD.${chrm}_${startpos}_${endpos}.recode.vcf
tabix tmp/gjfallsites_repRemoved.noAD.${chrm}_${startpos}_${endpos}.recode.vcf.gz

gjfs=tmp/gjfallsites_repRemoved.noAD.${chrm}_${startpos}_${endpos}.recode.vcf


#$poolAD and $gjf contains invar and vars that are at least with some levels of quality. Thus a ./. can be viewed as missing rather than an uncalled 0/0. 
bcftools merge $gjfs $poolAD -O z -o 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz
tabix 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz


python3 $parentDIR/useful_scripts/anc_state_infer.py  -I 116combineVCFs/116.${chrm}_${startpos}_${endpos}.vcf.gz -O 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf -Alist $TXTDIR/GJF.popline.txt

python3 $parentDIR/useful_scripts/lineage_polymorphisms.py -I 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf.gz -A $TXTDIR/RJF.popline.txt -Ao lineage_vcfs/RJF.polymorphisms.${chrm}_${startpos}_${endpos}.vcf -B $TXTDIR/DC.98.DL.popline.txt -Bo lineage_vcfs/DC.98.polymorphisms.${chrm}_${startpos}_${endpos}.vcf

#gunzip 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf.gz
#bgzip 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf
#tabix 116combineVCFs/swapped.${chrm}_${startpos}_${endpos}.vcf.gz

#parallel file=$TXTDIR/DC.98.DL.popline.txt or RJF.popline.txt 
vcftools --vcf lineage_vcfs/DC.98.polymorphisms.${chrm}_${startpos}_${endpos}.vcf \
        --indv INDV \
        --non-ref-ac 1 --recode \
        --out indv_vcfs/INDV.DC_poly.swapped.${chrm}_${startpos}_${endpos}
bgzip indv_vcfs/INDV.DC_poly.swapped.${chrm}_${startpos}_${endpos}.recode.vcf
tabix indv_vcfs/INDV.DC_poly.swapped.${chrm}_${startpos}_${endpos}.recode.vcf.gz







