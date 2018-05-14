import sys
import numpy as np
import stats_utilities as st
from joblib import Parallel, delayed
import multiprocessing
import glob

def get_eqtl(in_f):
        f_open = open(in_f)
        next(f_open)
        eqtl = {}
        for line in f_open:
                items = line.split()
                if items[1] not in eqtl:
                        eqtl[items[1]] = [[],[]]
                        eqtl[items[1]][0].append(line)
                        eqtl[items[1]][1].append(float(items[4]))
                        #print eqtl
                        #break
                else:
                        eqtl[items[1]][0].append(line)
                        eqtl[items[1]][1].append(float(items[4]))
        return eqtl


def BY_adj(eqtls):
        min_BY = []
        for gene in eqtls:
                pvalues = eqtls[gene][1]
                temp = st.pvalues(np.asarray(pvalues))
                BY_p = st.pvalues.BY(temp)
                eqtls[gene].append(BY_p)
                min_BY.append(min(BY_p))
        return eqtls,min_BY            

def process(in_f):
        print("getting eQTLs....")

        eqtls = get_eqtl(in_f)

        print("BY correction....")
        eqtls_BY, min_BY = BY_adj(eqtls)


        print("BH correction....")
        min_BY_p =  st.pvalues(np.asarray(min_BY))
        BY_BH_p = st.pvalues.BH(min_BY_p)

        max_p = max(i for i in BY_BH_p if i < 0.1)
        num = (BY_BH_p<0.1).sum()
        num1 = (BY_BH_p<0.05).sum()
        num2 = (BY_BH_p<0.01).sum()
        print(in_f + "\n number of eGenes under 0.1 level: %d" %num + "\n number of eGenes under 0.05 level: %d" %num1 + "\n number of eGenes under 0.01 level: %d" %num2)

        max_p_index = BY_BH_p.tolist().index(max_p)
        BY_cutoff = min_BY[max_p_index]
        BY_BH_p = BY_BH_p.tolist()
        out_name = in_f.split(".")[0]+"_BY_BH_sig_0.1.out"
        out_n = open(out_name,"w")
        out_n.write("SNP        gene    beta    t-stat  p-value FDR     BY_pvalue       BY_BH_value\n")
        write_sig(out_n,BY_cutoff,eqtls_BY,BY_BH_p)





def write_sig(out_n,BY_cutoff,eqtls,BY_BH_p):
        for t, gene in enumerate(eqtls):
                for n, i in enumerate(eqtls[gene][2]):
                        if i < BY_cutoff:
                                line = eqtls[gene][0][n].strip()+"\t"+str(eqtls[gene][2][n])+ "\t" + str(BY_BH_p[t]) +"\n"
                                out_n.write(line)
        out_n.close()

#files = sorted(glob.glob("./*covariate.out"))
#print(files)
#for i in files:
#        process(i)
#num_cores = multiprocessing.cpu_count()
#Parallel(n_jobs=num_cores)(delayed(process)(i) for i in files)


process(sys.argv[1])
