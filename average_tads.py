import argparse


import matplotlib
matplotlib.use('Agg') # runs headless

import numpy as np
import scipy as sp
import cooler
import klib

import matplotlib.pyplot as plt
import pdb
import sys
import os



parser=argparse.ArgumentParser(description='Create average TAD tip plot',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-ic',help='input matrix (cool format)',dest='incool',type=str,required=True)
parser.add_argument('-it',help='input TAD file (format: Chrom\tStart_pos\tEnd_pos',dest='intad',type=str,required=True)
parser.add_argument('-o',help='output prefix (if supplied, saves scaling plot and estimated coefficient)',dest='outprefix',type=str, required=True)

args=parser.parse_args()

tadfile = args.intad
coolfile = args.incool
outprefix = args.outprefix

def main():
  
    c=cooler.Cooler(coolfile)
    chroms=c.chromnames[:-1]

    winsize=25

    m_d=[]
    m_oe=[]
    for chrom in chroms:
        print(chrom,file=sys.stderr)
        ex = c.extent(chrom)
        d = c.matrix()[ex[0]:ex[1],ex[0]:ex[1]]
        n=d.shape[0]
        d[np.tril_indices(n,1)] = np.nan
        oe = obs_over_exp(d)

        with open(tadfile,'r') as fh:
            for x in fh:
                x = x.rstrip("\n").split("\t")
                if x[0]==chrom:
                    ex2 = c.extent(x[0]+':'+x[1]+'-'+x[2])
                    ex2 = (ex2[0]-ex[0],ex2[1]-ex[0])
                    assert ex2[0]<ex2[1]
                    if ex2[1]+winsize+1<=n and ex2[0]-winsize>=0:
                        d2 = d[ex2[0]-winsize:ex2[0]+winsize+1,ex2[1]-winsize:ex2[1]+winsize+1]
                        oe2 = oe[ex2[0]-winsize:ex2[0]+winsize+1,ex2[1]-winsize:ex2[1]+winsize+1]
                        m_d.append(d2)
                        m_oe.append(oe2)


    m_d=np.dstack(m_d)
    m_oe=np.dstack(m_oe)
    m_loe=np.log(m_oe)
    m_loe[~np.isfinite(m_loe)]=np.nan

    mean_d = np.nanmean(m_d,2)
    mean_oe = np.nanmean(m_oe,2)
    mean_loe = np.nanmean(m_loe,2)
    
    plt.figure()
    klib.heatmap(mean_d)
    np.save(outprefix+'_mean_d',mean_d)
    plt.savefig(outprefix+'_mean_d.svg')
    plt.savefig(outprefix+'_mean_d.png')

    plt.figure()
    klib.heatmap(mean_oe)
    np.save(outprefix+'_mean_oe',mean_oe)
    plt.savefig(outprefix+'_mean_oe.svg')
    plt.savefig(outprefix+'_mean_oe.png')    


    plt.figure()
    klib.heatmap(mean_loe,cmap='worb')
    np.save(outprefix+'_mean_loe',mean_loe)
    plt.savefig(outprefix+'_mean_loe.svg')
    plt.savefig(outprefix+'_mean_loe.png')





def obs_over_exp(d):
    n=d.shape[0]
    result=np.zeros((n,n))

    for i in range(2,500):
        print(i,file=sys.stderr)
        diag=np.diag(d,i)
        newdiag = diag/np.mean(diag[np.isfinite(diag)])
        result += np.diag(newdiag,i)

    result[np.tril_indices(n,1)]=np.nan
    result[np.triu_indices(n,500)]=np.nan

    return result


main()




