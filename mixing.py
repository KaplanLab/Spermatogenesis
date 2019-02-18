

import argparse

parser=argparse.ArgumentParser(description='Simulate mixing and evaluate GC strength',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-ps',help='PS cool file',dest='PS_path',type=str,required=True)
parser.add_argument('-s',help='S cool file',dest='S_path',type=str,required=True)
parser.add_argument('-o',help='output prefix',dest='out_prefix',type=str,required=True)
parser.add_argument('-r',help='resolution',dest='res',type=str,default='256000')
parser.add_argument('-p',help='purity (fraction)',dest='M_purity',type=float,default=0.9)
parser.add_argument('-d',help='use top and bottom 1/d of GC eigenvalues to calculate GC strength',dest='d',type=int,default=4)
parser.add_argument('-diag',help='only use these diagonals for GC strength calculation [FROM TO]. zero-based, inclusive. negative value x represent diagonal n+x. default is [1 n].',dest='use_diags',type=int,nargs=2)
parser.add_argument('-sgc',help='Use GC calculated from S for all datasets',dest='use_sgc',action='store_true')
parser.add_argument('-show',help='Show interactive plots',dest='show',action='store_true')
parser.add_argument('-seed',help='Randomization seed',dest='seed',type=int,default=1)

args=parser.parse_args()

PS_path = args.PS_path
S_path = args.S_path
out_prefix = args.out_prefix
res = args.res
M_purity = args.M_purity
d = args.d
use_diags = args.use_diags
use_sgc = args.use_sgc
show = args.show
seed = args.seed

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
import klib


def main():
    

    fh_stats = open(out_prefix+'_stats.tab','w')

    print ("PS_path:",PS_path,file=fh_stats)
    print ("S_path:",S_path,file=fh_stats)
    print ("resolution:",res,file=fh_stats)
    print ("purity:",M_purity,file=fh_stats)
    print ("d:",d,file=fh_stats)
    print ("use_diags:",use_diags,file=fh_stats)
    print ("use_sgc:",use_sgc,file=fh_stats)
    print ("seed:",seed,file=fh_stats)
    
    PS_cool = cooler.Cooler(PS_path+'::resolutions/'+res)
  
    x_extent = PS_cool.extent('chrX')
    n = x_extent[0]

    PS = PS_cool.matrix()[:n,:n] # PS will include all chromosomes up to chrX
    PS_x = PS_cool.matrix()[x_extent[0]:x_extent[1]+1,x_extent[0]:x_extent[1]+1] # PS_x is the matrix of chrX
    PS_triusum = np.nansum(np.triu(PS_cool.matrix(balance=False)[:n,:n],1))

    bin_chr=PS_cool.bins()['chrom'][:n].values.astype(str)

    # construct compartmentless ("expected") matrix based on PS

    PS_expected = klib.expected_matrix(PS,mode='cistrans',bin_chr=bin_chr)

    S_cool = cooler.Cooler(S_path+'::resolutions/'+res)
    S = S_cool.matrix()[:n,:n]
    S_x = S_cool.matrix()[x_extent[0]:x_extent[1]+1,x_extent[0]:x_extent[1]+1]
    S_triusum = np.nansum(np.triu(S_cool.matrix(balance=False)[:n,:n],1))

    
    prng = np.random.RandomState(seed)

    C = sample_matrix(S*(1-M_purity),int(PS_triusum*(1-M_purity)),prng) + sample_matrix(PS_expected*M_purity,int(PS_triusum*M_purity),prng) # We create a mixed matrix with a fraction 1-M sampled from S and fraction M sampled from PS

    S_expected = klib.expected_matrix(S,mode='cistrans',bin_chr=bin_chr)
    C_expected = klib.expected_matrix(C,mode='cistrans',bin_chr=bin_chr)

    PS_x_expected = klib.expected_matrix(PS_x,mode='cis',bin_chr=bin_chr)
    S_x_expected = klib.expected_matrix(S_x,mode='cis',bin_chr=bin_chr)
   
   
    PS_gc=klib.genomic_compartments(PS,PS_expected)
    S_gc=klib.genomic_compartments(S,S_expected)
    C_gc=klib.genomic_compartments(C,C_expected)
    PS_x_gc = klib.genomic_compartments(PS_x,PS_x_expected)
    S_x_gc = klib.genomic_compartments(S_x,S_x_expected)


    alt_gc = None
    if use_sgc:
        PS_gc = S_gc
        C_gc = S_gc
        PS_x_gc = S_x_gc

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(PS,PS_expected,PS_gc,klib.create_mask(bin_chr,'trans'),d=d)
    plot(PS,gc_str,PS_gc,bottomq,topq,LOE,sorted_LOE,'PS_trans',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(PS,PS_expected,PS_gc,klib.create_mask(bin_chr,'cis'),d=d)
    plot(PS,gc_str,PS_gc,bottomq,topq,LOE,sorted_LOE,'PS_cis',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(S,S_expected,S_gc,klib.create_mask(bin_chr,'trans'),d=d)
    plot(S,gc_str,S_gc,bottomq,topq,LOE,sorted_LOE,'S_trans',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(S,S_expected,S_gc,klib.create_mask(bin_chr,'cis'),d=d)
    plot(S,gc_str,S_gc,bottomq,topq,LOE,sorted_LOE,'S_cis',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(C,C_expected,S_gc,klib.create_mask(bin_chr,'trans'),d=d)
    plot(C,gc_str,S_gc,bottomq,topq,LOE,sorted_LOE,'mixed_trans',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(C,C_expected,S_gc,klib.create_mask(bin_chr,'cis'),d=d)
    plot(C,gc_str,S_gc,bottomq,topq,LOE,sorted_LOE,'mixed_cis',out_prefix,fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(PS_x,PS_x_expected,S_x_gc,d=d)
    print ('PS_chrX gc strength:',gc_str,file=fh_stats)

    gc_str, bottomq, topq, LOE, sorted_LOE = klib.gc_strength(S_x,S_x_expected,S_x_gc,d=d)
    print ('S_chrX gc strength:',gc_str,file=fh_stats)

    fh_stats.close()

    if show:
        
        plt.show()



def plot(mat,gc_str,gc,bottomq,topq,LOE,sorted_LOE,name,out_prefix,fh_stats):

    plt.figure(figsize=(7,10))
    ax=plt.subplot2grid((4,1),(0,0),rowspan=3)
    plt.title(name+' GC strength: %.2f' % gc_str)
    plt.imshow(np.log(mat))
    plt.subplot2grid((4,1),(3,0),sharex=ax)
    plt.plot(gc)
    plt.plot(topq,gc[topq],'.')
    plt.plot(bottomq,gc[bottomq],'.')
    plt.savefig(out_prefix+'_'+name+'_GCstr1.png')

    plt.figure(figsize=(15,5))
    plt.subplot(1,3,1)
    plt.imshow(LOE,vmin=-2,vmax=2,cmap='RdBu')
    plt.title('log Obs/Exp '+name)
    plt.subplot(1,3,2)
    plt.imshow(sorted_LOE,vmin=-2,vmax=2,cmap='RdBu')
    plt.title('log Obs/Exp '+name+' sorted by GC')
    q1 = len(topq)
    q2 = sorted_LOE.shape[0]-q1
    draw_square(0,0,q1,color='blue')
    draw_square(q2,q2,q1,color='blue')
    draw_square(0,q2,q1,color='purple')
    draw_square(q2,0,q1,color='purple')
    plt.xlim(0,sorted_LOE.shape[0])
    plt.ylim(sorted_LOE.shape[0],0)
    plt.colorbar()
    plt.subplot(1,3,3)
    plt.savefig(out_prefix+'_'+name+'_GCstr2.png')

    print (name,'gc strength:',gc_str,file=fh_stats)


def draw_square(top_left_x, top_left_y, length, lw=2, color='orange'):
    plt.plot([top_left_x,top_left_x + length],[top_left_y,top_left_y], lw=lw,color=color)
    plt.plot([top_left_x,top_left_x + length],[top_left_y + length,top_left_y+length], lw=lw,color=color)
    plt.plot([top_left_x,top_left_x],[top_left_y,top_left_y + length], lw=lw,color=color)
    plt.plot([top_left_x + length,top_left_x + length],[top_left_y, top_left_y + length],lw=lw,color=color)



def sample_matrix(A,readnum,prng):
    '''
    Converts A into probability matrtix and samples readnum reads multinomially according to the probabilities of A. Symmetry is kept. The resulting matrix is renormalized so its sum is equal that of the original A (important because A may be balanced so might have a weird sum).
    '''

    Asum = np.nansum(A)
    A=np.nan_to_num(np.triu(A,0))    
    Av=A.view().reshape((-1,))
    Av[:] = prng.multinomial(readnum,Av/np.sum(Av))
    A = A + np.triu(A,1).T
    A = (Asum/np.sum(A))*A

    return A





main()

