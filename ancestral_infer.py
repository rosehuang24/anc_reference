#this script aims to infer the ancestral allele based on the outgroup genotyp information. However, this script is highly specialized for my project and should be carefully ediited when applied to another.
#The input files were obtained by GATK4.1.2 with --all-site parameter when calling vcfs. Here we did not joint call but rather merge data with bcftools. And due to the presence of invariants, we have to abandon AD information, resulting in an incomplete vcfs (compare to convetional vcfs.)
#prior ro merging inndividual vcfs, we also filtered them with primitive standards, such as quality (mean depth, minQual)and max missing genotype calling (2 indvs for GJF.)
#Also, the merged vcfs were also processed by vcf to remove repetitve regions and to contain at least one alternate allele.
#Therefore I don't have to consider the invariants sites possessed by all samples (which is totally a waste of calculation.)


import argparse
import gzip
parser = argparse.ArgumentParser()
parser.add_argument("-I","--input", help="input file, zipped or unzipped vcf",required=True)
parser.add_argument("-O","--output", help="output vcf file name, can be ended with .gz",required=True)
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


if args.output.endswith(".gz"):
    outh = gzip.open(args.output, 'wt')
else:
    outh = open(args.output, 'w')

indices_pop=[]
#swap_dict={"1/1":"0/0","0/0":"1/1","0/1":"0/1","./.":"./.",".":"./."}

with (gzip.open if args.input.endswith(".gz") else open)(args.input, "rt") as inh:
    for lines in inh:
        if lines.startswith("##"):
            outh.write(lines)
        elif lines.startswith("#CHROM"):
            outh.write(lines)
            line=lines.strip().split()
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
            miss=int(ls.count("./."))-int(anc_ls.count(./.))
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
                    outh.write("\t".join(line[:3])+"\t"+line[4]+"\t"+line[3]+"\t"+"\t".join(line[5:7])+"\tStatsSwapped\tGT\t"+"\t".join(ls))
                elif anc_ls.count("0/0")>=7 and anc_ls.count("1/1") == 0:
                    outh.write(lines)
            #print(anc_ls.count("0/0"))




inh.close()
outh.close()
