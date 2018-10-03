from collections import OrderedDict, Counter
import sys
import glob
import re
from os.path import basename
import os
def process(bed):
        '''
        input: bed file from RFMix output
        output: [{region: pop, region: pop}, ...22 dicts in total]
        '''
        bed_a = open(bed)
        count = 0
        tracts = []
        chr_dict = {}
        for line in bed_a:
                items = line.split()
                if count != int(items[0]):
             
                        if bool(chr_dict):
                                tracts.append(chr_dict)        
                                chr_dict = {}
                                count = int(items[0])
                                chr_dict[tuple([int(items[1]), int(items[2])])] = items[3]
                        
                        else:
                                chr_dict[tuple([int(items[1]), int(items[2])])] = items[3]
                                count = int(items[0])
                else:
                        chr_dict[tuple([int(items[1]), int(items[2])])] = items[3]
                
        tracts.append(chr_dict)
        tracts_order = [OrderedDict(sorted(dic.items(),key=lambda t: t[0][0])) for dic in tracts]
        return(tracts_order)

               

def get_per_gene(tracts_order, gene_interval, gene_chr):
        '''
        input: tracts for sampleX, gene interval, gene chr
        output:
        '''

        tracts_dict = tracts_order[int(gene_chr)-1]
        percent = 0
        function_interval = 0
        for interval in tracts_dict.keys():
                ancestry = 1 if tracts_dict[interval] == "AFR" else 0
                if gene_interval[1] <= interval[0]:
                        break

                elif gene_interval[1] <= interval[1]:
                        if gene_interval[0] <= interval[0]:
                                temp = gene_interval[1] - interval[0]
                                percent = percent + temp * ancestry
                                function_interval = function_interval + temp                              # break
                        elif gene_interval[0] >= interval[0]:
                                temp = gene_interval[1] - gene_interval[0]
                                function_interval = function_interval + temp
                                percent = percent + temp * ancestry
                        break
                elif gene_interval[1] >= interval[1]:
                        if gene_interval[0] > interval[1]:
                                continue
                        else:
                                if gene_interval[0] >= interval[0]:
                                        temp = interval[1] - gene_interval[0] 
                                        percent = percent + temp * ancestry
                                        function_interval = function_interval + temp
                                elif gene_interval[0] < interval[0]:
                                        temp = interval[1] - interval[0] 
                                        function_interval = function_interval + temp
                                        percent = percent + temp * ancestry
                                continue
                                       
        if function_interval == 0:
                res = "NA"
        else:
                res = float(percent)/float(function_interval)
               
               
        return res

def get_pos(bed, tracts_order_a, tracts_order_b, genesPer):
        chrs = [str(i) for i in range(1, 23)] 
        '''
        input: bed file containing gene positon, tract1 for sampleX, tract2 for sampleX, genepercentage dict
        output: genepercentage dict after adding the percentage for this sampleX
        '''
        f = open(bed)
        for line in f:

                if line[0] == "#":
                        continue 
                items = line.split()
                chrGene = items[0]
                if chrGene in chrs:
                        start = int(items[1])
                        end = int(items[2])
                        percent1 = get_per_gene(tracts_order_a,[start, end], chrGene)
                        percent2 = get_per_gene(tracts_order_b,[start, end], chrGene)
                        if "NA" in [percent1, percent2]:
                                percent_gene = "NA"
                        else:                                             
                                percent_gene = (percent1+percent2)/2
                        if tuple(items[:4]) in genesPer:
                                genesPer[tuple(items[:4])].append(percent_gene)
                        else:
                                genesPer[tuple(items[:4])] = [percent_gene]
                        
                else:
                        continue

        return genesPer



def main(args):
        bed = args.bed
        working_dir = args.dirs
        beds = glob.glob(working_dir + "*A.bed")
        samples = ["chr", "start", "end", "gene"]
        genePer = {}

        for rfbed in beds:
                       
                print(rfbed)
                tracts_order_a = process(rfbed)
                samples.append(re.sub(r'_pipe_A.bed', '', basename(rfbed)))
                tracts_order_b = process(re.sub(r'A.bed', 'B.bed', rfbed))
                genesPer = get_pos(bed, tracts_order_a, tracts_order_b, genePer)
                

        res = open(working_dir + "AFR_percentage.txt","w")
        res.write("\t".join(samples)+"\n")
        for i in genesPer.keys():
                item = list(i) +  [str(s) for s in genesPer[i]]
                
                res.write("\t".join(item) + "\n")


        os.system("(head -n 2 " + working_dir + "AFR_percentage.txt && tail -n +3 " + working_dir + "AFR_percentage.txt | sort -k1,1 -k2,2n) > " + working_dir + "AFR_percentage.sorted.txt ")

if __name__ == '__main__':
        parser = argparse.ArgumentParser()
        parser.add_argument('--bed',help='bed file with gene postion')
        parser.add_argument('--dirs',help='working directory')
        args = parser.parse_args()
        main(args)
