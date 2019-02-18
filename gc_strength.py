

import argparse

parser=argparse.ArgumentParser(description='Calculate GC strength with respect to a given genomic compartment vector',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i',help='input cool file',dest='infile',type=str,required=True)
parser.add_argument('-gc',help='optional secondary cool file used to calculate genomic compartment vector',dest='gcfile',type=str)
parser.add_argument('-c',help='calculate only cis gc_strength for this chromosome',dest='singlechr',type=str)
parser.add_argument('-r',help='resolution',dest='res',type=str,default='250000')
parser.add_argument('-d',help='use top and bottom 1/d of GC eigenvalues to calculate GC strength',dest='d',type=int,default=4)


args=parser.parse_args()

infile=args.infile
gcfile=args.gcfile
singlechr=args.singlechr
res=args.res
d=args.d


import numpy as np
import scipy as sp
import cooler
import sys
import klib
import matplotlib.pyplot as plt
import pdb


def main():
    
    incool = cooler.Cooler(infile+'::resolutions/'+res)
    
    if not singlechr:

        A = incool.matrix()[:]
        bin_chr=incool.bins()['chrom'][:].values.astype(str)

    else:
        c_extent = incool.extent(singlechr)
        A = incool.matrix()[c_extent[0]:c_extent[1],c_extent[0]:c_extent[1]]
        bin_chr=incool.bins()['chrom'][c_extent[0]:c_extent[1]].values.astype(str)


    A_expected = klib.expected_matrix(A,mode='cistrans',bin_chr=bin_chr)

    if not gcfile:
        gc=klib.genomic_compartments(A,A_expected)
    else:
        scool=cooler.Cooler(gcfile+'::resolutions/'+res)
        
        if not singlechr:
            S = scool.matrix()[:]
        else:
            S = scool.matrix()[c_extent[0]:c_extent[1],c_extent[0]:c_extent[1]]

        gc=klib.genomic_compartments(S,klib.expected_matrix(S,mode='cistrans',bin_chr=bin_chr))
  

    gc_str = klib.gc_strength(A,A_expected,gc,klib.create_mask(bin_chr,'trans'),d=d)[0]
    print ("trans GC strength:",gc_str)
    gc_str = klib.gc_strength(A,A_expected,gc,klib.create_mask(bin_chr,'cis'),d=d)[0]
    print ("cis GC strength:",gc_str)







main()

