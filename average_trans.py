

import argparse

parser=argparse.ArgumentParser(description='Calculate averaged trans (interchromosomal) interaction map',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-ps',help='input cool file',dest='PS_path',type=str,required=True)
parser.add_argument('-o',help='output prefix',dest='out_prefix',type=str,required=True)
parser.add_argument('-r',help='resolution',dest='res',type=str,default='256000')

args=parser.parse_args()

PS_path = args.PS_path
out_prefix = args.out_prefix
res = args.res


import matplotlib
import numpy as np
import scipy as sp
import cooler
import scipy.spatial
import matplotlib.pyplot as plt
import sklearn.decomposition
import sys
import pdb
import time
import PIL
    

def main():
    
    m=500

    PS_cool = cooler.Cooler(PS_path+'::resolutions/'+res)
    n = PS_cool.extent('chrX')[0]
    bin_chr=PS_cool.bins()['chrom'][:n].values.astype(str)
    chrs=np.unique(bin_chr)
    c=len(chrs)

    PS = PS_cool.matrix()[:n,:n]
    print (np.nansum(PS))

    M=np.zeros((m,m,int(c*(c-1))))
    M[:]=np.nan
    i=0
    for c1 in chrs:
        for c2 in chrs:
            if c1 != c2:
                print(c1,c2)
                D = PS[bin_chr==c1,:][:,bin_chr==c2]
                x = np.array(PIL.Image.fromarray(D).resize((m,m)))
                M[:,:,i] = np.clip(x,None,np.nanpercentile(x,99.9))
                i+=1

    r=np.nanmean(M,2)
    
    r = np.log(r/np.nanmean(r))
    plt.figure(figsize=(10,8))


    plt.imshow(r,vmin=-0.4,vmax=1,cmap="seismic",aspect='equal')
    plt.xlim(0,m)
    plt.ylim(m,0)
    plt.colorbar()
    plt.savefig(out_prefix+'.png')
    plt.savefig(out_prefix+'.svg')


main()
