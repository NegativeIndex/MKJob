#!/usr/bin/env python3
import os
import re
import subprocess
import glob
import numpy as np
import shutil 
import sys
sys.path.insert(0,'/Users/wdai11/function')
import my_output as my
import flatten
import itertools

######################### 
# global variables
##########################
class common:
    fname="rad.py"     
    otherfiles=()
    core=32
    memory=50
    queue='all.q'
   
######################### 
# generate one job
##########################
def find_line_numbers(lines,regexp):
    numbers=[]
    for i,line in enumerate(lines):
        # print(line)
        if re.match(regexp, line):
            numbers.append(i)
    return numbers

######################### 
# generate one job
##########################
def other_steps(dname):
    
    # copy files
    for ff in common.otherfiles:
        print('Copy file '+ff)
        shutil.copy(ff,dname)

    # generate job file
    os.chdir(dname)
    fname=common.fname
    subprocess.call(["mkjob-mpi",fname])
    jobfiles=glob.glob("dwt*job")
    ss_core="-c {}".format(common.core)
    ss_mem="-m {}".format(common.memory)
    ss_q="-q {}".format(common.queue)
    subprocess.call(["modify-meep-jobfile.py"]+jobfiles
                    +[ ss_core,ss_mem,ss_q  ])       
    # generate job.begin
    open("job.begin", 'a').close()
    print("Generate job.begin")

######################### 
# check two files
##########################
def check_two_files():
    with open(common.fname) as f:
        file1=f.read().splitlines()

    with open("create_3D_jobs.py") as f:
        file2=f.read().splitlines()

    file1=[line.lstrip().rstrip() for line in file1 
           if re.search('=',line)]
    file2=[line.lstrip().rstrip()
           for line in file2 if re.search('=',line)]
    
    sameline=[line for line in file1 if line in file2]

    print('-----------------------')
    for line in sameline:
        print(line)
    print('-----------------------')
    userInput = input('The two files have such commone lines. y/n?  ');

    if userInput!='y':
        sys.exit('Something is wrong. Modify the file')
    
########################## 
# main function
##########################
# collect job file parameters
check_two_files()

# prepare py file
path=os.path.abspath("./")
fname=common.fname

lines=[]
with open(fname, "rt") as fin:
    fin=open(fname, 'rt')
    lines=fin.readlines()

line1=find_line_numbers(lines,'class common')
line2=find_line_numbers(lines,r'^##########################')

lmin=line1[0]
i=0
while line2[i]<lmin: i+=1
lmax=line2[i]
############################
# all the parameters
pols=["Ex",]
geoms=["Empty","Wg"]
fcens=np.arange(0.25,0.26,0.02)
dfs=np.ones(fcens.size)*0.02
freqs=list(zip(fcens,dfs))

gparas=[(nDBR,a,t,dx,dy,dz)
        for nDBR in np.array([1,3,5,7,10])
        for a in np.array([1.0,])
        for t in np.arange(0.8,0.9,1.0)
        for dx in np.array([0.0,])
        for dy in np.array([0.0,])
        for dz in np.arange(0.21,0.30,0.02)  ]
       
paras1=list(itertools.product(("Ex",), ("Empty",), freqs, 
                              ( (0,)*len(gparas[0]), )   ))
paras1=[flatten.flattenArrayN(elem) for elem in paras1 ]

paras2=list(itertools.product(pols, geoms[1:2], freqs, gparas ))
paras2=[ flatten.flattenArrayN(elem) for elem in paras2 ]

paras=paras1+paras2
# paras=paras1

# [print(elem) for elem in paras]

######################
njobs=len(paras)
print('{} jobs are generated'.format(njobs))
print('Queue {}, {} cores, {}G memery'.format(
    common.queue, common.core,common.memory))
userInput = input('Any key please, Ctrl-c to quit')

count=0
for para in paras:
    os.chdir(path)
    
    pol,geom,fcen,df,nDBR,a,t,dx,dy,dz=para

    sig='{}_F{:0.3f}'.format(geom,fcen)
    if geom != 'Empty':
        sig='{}_{}_F{:0.3f}'.format(geom,pol,fcen)
        pstring='_L{:d}_A{:0.2f}T{:0.2f}'.format(nDBR,a,t)
        dstring='_DX{:0.2f}Y{:0.2f}Z{:0.2f}'.format(dx,dy,dz)
        sig+=pstring+dstring
    
    dname=sig
    print(sig)
    if not os.path.exists(dname):
        os.makedirs(dname)
        # write the new py file
        fname2=os.path.join(dname,fname)
        for i in range(lmin,lmax):
            line=lines[i];
            # print(line)
            line=re.sub(r'^(\s+geom=).*$',r'\g<1>"{:s}"'.format(geom), line)
            line=re.sub(r'^(\s+pol=).*$',r'\g<1>"{:s}"'.format(pol), line)
            line=re.sub(r'^(\s+fcen=).*$', r'\g<1>{:0.3f}'.format(fcen), line)
            line=re.sub(r'^(\s+df=).*$', r'\g<1>{:0.3f}'.format(df), line)
            line=re.sub(r'^(\s+nDBR=).*$', r'\g<1>{:d}'.format(int(nDBR)), line)
            line=re.sub(r'^(\s+a=).*$', r'\g<1>{:0.3f}'.format(a), line)
            line=re.sub(r'^(\s+t=).*$', r'\g<1>{:0.3f}'.format(t), line)
            line=re.sub(r'^(\s+dx=).*$', r'\g<1>{:0.3f}'.format(dx), line)
            line=re.sub(r'^(\s+dy=).*$', r'\g<1>{:0.3f}'.format(dy), line)
            line=re.sub(r'^(\s+dz=).*$', r'\g<1>{:0.3f}'.format(dz), line)
            # print(line)
            lines[i]=line

        fout=open(fname2, 'wt')
        fout.write("".join(lines))
        fout.close()
        print("Generate {:s}".format(fname2))

        other_steps(dname)
        count+=1
        print('Generated {} jobs, {} to do'.format(count, njobs-count))
        print('-------------------')

print('{} jobs are generated'.format(count))

