import os
import sys

sds = []
count = 0

# For each job batch submitted
for i in range(0,800):

    os.system("mv ../../step2/wd/job"+str(i)+"/slurm* ../wd/job"+str(i)+"/slurm.out")

    filename = "../wd/job"+str(i)+"/slurm.out"
    f  = open(filename)
    A  = f.readlines()
    sz = len(A)
    f.close()
    for j in range(sz):
        if "excitations" in A[j]:

            sds.append( [] )

            b = A[j].split("],")   
            sds[count].append( [  b[0].strip().split("[")[2] ] ) 
            for k in range(1,len(b)):
                c = b[k].strip().split("[")
                g = c[1].strip().split("]")
                sds[count].append([g[0]])
            count += 1          

unique_sds = []
for i in range(len(sds)):
    for j in range(len(sds[i])):
        if sds[i][j] not in unique_sds:
            unique_sds.append(sds[i][j])    
print(unique_sds)
