
import argparse
import numpy as np
import pdb
import sys
import klib
import matplotlib
import matplotlib.pyplot as plt
import cooler

parser=argparse.ArgumentParser(description='visualize point interactions on Hi-C matrix',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# usage example: python -ic PS_pooled.mcool -il PS_loops.bedpe -g chr2:48e6-55e6 -r 20000

parser.add_argument('-ic',help='input matrix (cool format)',dest='incool',type=str,required=True)
parser.add_argument('-il',help='input point interactions (bedpe format (Vals not used): Chrom1\tStart_pos1\tEnd_pos1\tVal1\tChrom2\tStart_pos2\tEnd_pos2\tVal2) ',dest='inloop',type=str,required=True)
parser.add_argument('-g',help='genomic location (e.g. chr2:48000000-55000000)',dest='genloc',type=str, required=True)
parser.add_argument('-r',help='Hi-C resolution',dest='res',type=int, required=True)

args=parser.parse_args()


def main():

    cool_file = args.incool
    loop_file = args.inloop
    vrange = args.genloc
    res = args.res
    
    vchr = vrange.split(":")[0]
    vstart,vend = vrange.split(":")[1].split("-")
    vstart,vend = float(vstart),float(vend)

    assert vstart%res==0
    assert vend%res==0
    n = int((vend-vstart)/res)
    D=np.zeros((n,n))


    print('loading loops...',file=sys.stderr)
    loops=[]
    with open(loop_file,'r') as fh:
        next(fh)
        for x in fh:
            c1,s1,e1,_,c2,s2,e2,_ = x.rstrip("\n").split("\t")

            if c1==vchr and int(s1)>=vstart and int(e1)<=vend and c2==vchr and int(s2)>=vstart and int(e2)<=vend:
                loops += [(s1,e1,s2,e2)]

    print('done',file=sys.stderr)


    cool_file = cool_file+'::/resolutions/'+str(res)
    c = cooler.Cooler(cool_file)
    ex = c.extent(vchr+':'+str(int(vstart))+'-'+str(int(vend)))
    D=c.matrix()[ex[0]:ex[1],ex[0]:ex[1]]

    

    plt.figure()
    klib.heatmap(D,log=True,clip=0.02)
    
    for loop in loops:
        s1,e1,s2,e2 = map(float,loop)
        x = int((s1-vstart-1)/res)-0.5
        y = int((s2-vstart-1)/res)-0.5
        w = int((e1-s1-1)/res)+1
        h = int((e2-s2-1)/res)+1
       
        plt.gca().add_patch(matplotlib.patches.Rectangle((y,x),h,w,linewidth=2,edgecolor='b',facecolor='none'))
    plt.show()


main()
