#!/bin/bash
for filename in /Volumes/fsmresfiles/Basic_Sciences/Pharm/Perera_Lab/Yizhen/v7-GTEx/GTEx_Analysis_v7_eQTL/*gene_pairs.txt.gz; do
	echo $filename;
	zgrep -c $ $filename;
	awk 'NR==FNR{c[$2]++;next}$1 in c{print $0}' GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_maf0.05_SM.bim <(gzip -dc $filename) >temp.txt;
	wc -l temp.txt;
#awk 'NR==FNR{c[$1]++;next} !($1 in c)' *paste reserve.snp.paste.convert.add > reserve.snp.paste.not.tested
#awk '{print $1}' reserve.snp.paste.not.tested | sort | uniq | wc
#awk -F'[._\t]' '{print $6 "." $2 "\t" $0 }' Liver.v7.signif_variant_gene_pairs.txt >Liver.v7.signif_variant_gene_pairs.paste.txt
	awk -F'[._\t]' '{print $6 "." $2 "\t" $0 }' temp.txt >temp.paste.txt;
	awk 'NR==FNR{c[$1]++;next}$2 in c{print $1}' temp.paste.txt reserve.snp.paste.convert.add  | wc -l;
done
