import os
import glob
import sys
import subprocess

from joblib import Parallel, delayed
import multiprocessing


condit = sys.argv[1]
files = glob.glob("../data/genotype/eQTL/condition1/*condition"+condit+"*covariate.txt")

filesLA = glob.glob("../data/genotype/eQTL/condition1/*condition"+condit+"*covariate_nopc.txt")

script = "#!/bin/bash\n \
#MSUB -A b1042\n \
#MSUB -l walltime=7:00:00\n \
#MSUB -l nodes=1:ppn=1\n \
#MSUB -j oe\n  \
#MSUB -q genomics\n \
#SUB -N run_peer\n \
module load R \n \
cd /projects/b1047/zhong/Hepatocyte_project/script\n"

def runEach(index):

        fn = files[index]
        prefix = fn.split("/")[-1].split(".")[0]
        scriptR = "Rscript ../../software/ancestry_pipeline-master/eqtl.matrixEL.r \
        ../data/genotype/eQTL/condition" + condit +"_expression_tmm_normalization_exp_v7_dosage_v7_2.txt \
        ../data/genotype/eQTL/hepatocyte_snp.position ../data/genotype/eQTL/condition" + condit + "_expression_tmm_normalization_exp_v7.txt \
        ../data/genotype/eQTL/Hepatocyte_gene_pos.txt "+ fn + " ../results/eQTL/condition1/" + prefix +"_dosage.out"

        scriptF =  script + scriptR

        out =  "submit"+str(index+10)+"_"+condit+".sh"
        print(out)
        FileO =  open(out,"w")
        FileO.write(scriptF)
        FileO.close()

        print("chmod +x ./" + out)
        subprocess.call("chmod +x ./" + out ,shell=True)
        print("msub " + out)
        subprocess.call("msub " + out,shell=True)


def runEachLA(index):
        fn = filesLA[index]
        prefix = fn.split("/")[-1].split(".")[0]
        scriptR = "Rscript ../../software/ancestry_pipeline-master/matrixEL_local_cis.r \
        ../data/genotype/eQTL/condition" + condit +"_expression_tmm_normalization_exp_v7_dosage_v7_2.txt \
        ../data/genotype/eQTL/hepatocyte_snp.position ../data/genotype/eQTL/condition" + condit + "_expression_tmm_normalization_exp_v7.txt \
        ../data/genotype/eQTL/Hepatocyte_gene_pos.txt "+ fn + \
        " ../results/eQTL/condition1/" + prefix +"_LA_dosage.out \
        ../data/genotype/eQTL/condition" + condit + "_expression_tmm_normalization_ancestry_v7.txt"
        scriptF =  script + scriptR
        out =  "submitLA"+str(index+10)+"_"+condit+".sh"
        print(out)
        FileO =  open(out,"w")
        FileO.write(scriptF)
        FileO.close()
        print("chmod +x ./" + out)
        subprocess.call("chmod +x ./" + out ,shell=True)
        print("msub " + out)
        subprocess.call("msub " + out,shell=True)


num_cores = multiprocessing.cpu_count()
Parallel(n_jobs=num_cores)(delayed(runEach)(i) for i in range(len(files)))
print(files)

Parallel(n_jobs=num_cores)(delayed(runEachLA)(i) for i in range(len(filesLA)))

