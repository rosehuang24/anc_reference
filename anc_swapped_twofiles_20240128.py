#===========================================================================================================================================================================
#2024/01/28
# This script was directly editted from the anc_infer.py with modification of output files
# The reasons behind this file is that the SNPeff annotation has a pre-defined amino acid sequence 
# and when the reference and alternate state alleles were swapped, and if the change is non-synonymous,
# the resultant AA change will match the pre-defined sequence and missclassified this mutaton as synonymous.
# Therefore, I modified this script so the output files are (with the satisfaction that at least 7 GJF are all homozygous for 0 or 1):
# 1) all sites of all sameples, if allele states are wrong, the GT were swapped, but not the REF/ALT part. (easy for SNPeff annotation and genload calculation)
# 2) only sites with wrong allele states. The REF/ALT were swapped, no GT or any individual.

##============================================================================================================================================================================

#this script aims to infer the ancestral allele based on the outgroup genotyp information. However, this script is highly specialized for my project and should be carefully ediited when applied to another.
#The input files were obtained by GATK4.1.2 with --all-site parameter when calling vcfs. Here we did not joint call but rather merge data with bcftools. And due to the presence of invariants, we have to abandon AD information, resulting in an incomplete vcfs (compare to convetional vcfs.)
#prior ro merging inndividual vcfs, we also filtered them with primitive standards, such as quality (mean depth, minQual)and max missing genotype calling (2 indvs for GJF.)
#Also, the merged vcfs were also processed by vcf to remove repetitve regions and to contain at least one alternate allele.
#Therefore I don't have to consider the invariants sites possessed by all samples (which is totally a waste of calculation.)

##================================================================================================================================================

import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("-I","--input", help="input file, zipped or unzipped vcf",required=True)
parser.add_argument("-aO","--alloutput", help="output vcf file name for all sites. Individual genotypes included and swapped if true. Allele states were unchanged. Can be ended with .gz",required=True)
parser.add_argument("-sO","--swappedoutput", help="output vcf file name for sites with wrong allele states. No individual genotypes information.  Allele states were changed. can be ended with .gz",required=True)
parser.add_argument("-A","--ancindv",  nargs='+', help="individual name(s) of ancestral group, seperated with space.",required=False)
parser.add_argument("-Alist","--anclist",  help="a file with individual name(s) of ancestral group, one per line",required=False)

outgroup_indv=[]
args = parser.parse_args()
if args.ancindv and args.anclist:
    print("\nError: Only one kind of ancestral information accepted.\n")
else:
    if args.ancindv:
        outgroup_indv=args.ancindv
    else:
        if args.anclist:
            for lines in open(args.anclist,'r'):
                outgroup_indv.append(lines.strip())
        else:
            print("\nError: Ancestral information required!\n")


if args.alloutput.endswith(".gz"):
    outah = gzip.open(args.alloutput, 'wt')
else:
    outah = open(args.alloutput, 'w')

if args.swappedoutput.endswith(".gz"):
    outsh = gzip.open(args.swappedoutput, 'wt')
else:
    outsh = open(args.swappedoutput, 'w')

indices_pop=[]
#swap_dict={"1/1":"0/0","0/0":"1/1","0/1":"0/1","./.":"./.",".":"./."}

with (gzip.open if args.input.endswith(".gz") else open)(args.input, "rt") as inh:
    for lines in inh:
        if lines.startswith("##"):
            outah.write(lines)
            outsh.write(lines)
        elif lines.startswith("#CHROM"):
            outah.write("INFO=<ID=StatsUnswapped,Number=0,Type=Flag,Description=\"Allele state change\""+"\n")
            outsh.write("INFO=<ID=StatsUnswapped,Number=0,Type=Flag,Description=\"Allele state change\""+"\n")
            outah.write(lines)
            line=lines.strip().split()
            outsh.write('\t'.join(line[0:9])+"\n")
            indices_pop = [int(line.index(x))-9 for x in outgroup_indv]

            '''outh.write('\t'.join(line[0:9]))
            for i in indices_pop:
                outh.write('\t'+line[i])
            outh.write('\n')'''
        elif not lines.startswith("#"):
            ls=[]
            anc_ls=[]

            for field in lines.strip().replace("|", "/").split()[9:]:

                gt=str(field).split(":")[0]
                ls.append(gt)#so the genotype in ls has the corresponding indices with the individual in "indices_pop" list.
            #for i in indices_pop:
            #   anc_ls.append(ls[int(i)])
            anc_ls=[ls[int(i)] for i in indices_pop]
            #print(anc_ls)
            miss=int(ls.count("./."))-int(anc_ls.count("./."))
            if "0/1" not in anc_ls and miss<12:

                if anc_ls.count("1/1")>=7 and anc_ls.count("0/0") == 0:
                    line=lines.strip().split()
                    for x in range(len(ls)):
                        if ls[x] == '1/1':
                            ls[x] = '0/0'
                        elif ls[x] == '0/0':
                            ls[x] = '1/1'
                        elif ls[x] == '.':
                            ls[x]='./.'
                #    print(ls)
                    outsh.write("\t".join(line[:3])+"\t"+line[4]+"\t"+line[3]+"\t"+"\t".join(line[5:7])+"\tStatsSwapped\tGT"+"\n")
                    outah.write('\t'.join(line[0:7])+"\tStatsUnswapped\tGT\t"+"\t".join(ls)+"\n")
                elif anc_ls.count("0/0")>=7 and anc_ls.count("1/1") == 0:
                    outah.write(lines)
            #print(anc_ls.count("0/0"))




inh.close()
outah.close()
outsh.close()
