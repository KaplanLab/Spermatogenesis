
import sys
import numpy as np
import scipy as sp

import sklearn.decomposition

import matplotlib
import matplotlib.pyplot as plt
import pdb

worb_cmap=matplotlib.colors.LinearSegmentedColormap.from_list('worb',colors=['white','orange','red',[0.5,0,0],'black'])
worb_cmap.set_bad([0.82,0.82,0.82])
plt.register_cmap(cmap=worb_cmap)


def heatmap(A,cmap="worb",clip=0,top_half=False,log=False,colorbar=True,chrs=None,vmin=None,vmax=None,imshow_interpolation="none",coords=False):
    '''
    clip_top clips top fraction of data
    top_half shows only top_half of matrix
    '''
    
    if clip>0:
        A=np.clip(A,-np.inf,np.percentile(A[~np.isnan(A)],(1.0-clip)*100))
    if log:
        A=np.log(1+A)
    
    plt.imshow(A,interpolation=imshow_interpolation,cmap=cmap,vmin=vmin,vmax=vmax)
    if chrs != None:
            
        for i in range(1,len(chrs)):
            if chrs[i]!=chrs[i-1]:
                plt.axvline(i-0.5,color="black")
                plt.axhline(i-0.5,color="black")

    if colorbar:
        plt.colorbar()



def genomic_compartments(A,expected_A):
    '''
    calculates genomic compartments by applying PCA to the correlation matrix of normalizing A/expected_A and taking the first eigenvector
    '''

    n=A.shape[0]

    OE = (A + 1e-7) / (expected_A + 1e-7) # add pseudocounts to solve zeros

    OE = np.clip(OE,None,np.nanpercentile(OE,99.9))

    valid = ~np.all(np.isnan(OE),0)
    clean_OE = OE[valid,:][:,valid]

    corrmat = np.corrcoef(clean_OE)

    pca=sklearn.decomposition.PCA(n_components=1)
    pca.fit(corrmat)
    comp=np.empty(n)
    comp[:]=np.nan
    comp[valid]=pca.components_[0]

    return comp


def expected_matrix(A,mode='cis',bin_chr=None,func=np.nanmedian):
    '''
    calculate "expected" matrix for Hi-C matrix A.
    in cis, this means setting the values of diagonal d to func(d).
    in trans, this means setting all values to func(all trans values).
    
    if mode='cis', assumes A is only a cis matrix
    if mode='cistrans', assumes A contains both cis and trans
    
    bin_chr: array giving the chr to which each bin belongs
    '''

    n=A.shape[0]

    if mode=='cis':
        B=np.zeros((n,n))

        for i in range(n):
            diag_ind=(range(i,n),range(0,n-i))
            x = func(A[diag_ind])

            B[diag_ind] = x
            B[diag_ind[::-1]] = x


    if mode=='cistrans':

        B=np.zeros((n,n))
        trans_mask = create_mask(bin_chr,'trans')
        B[trans_mask]=func(A[trans_mask])                                                                                                                               
        # each chromosome is normalized separately
        chrs=np.unique(bin_chr)
        for c in chrs:
            b = np.where(bin_chr==c)[0]
            m=b.shape[0]
            for i in range(m):
                inds = (b[i:m],b[0:m-i])
                x = func(A[inds])
                B[inds] = x
                B[inds[::-1]] = x

    return B





def create_mask(bin_chr,mode='cis',diags=None):
    '''
    mode 'cis': return a matrix where all cis bins are True and all trans are false
    mode 'trans': return a matrix where all trans bins are True and all cis are false
    
    bin_chr: array giving the chr to which each bin belongs
    diags
    '''

    n=bin_chr.shape[0]
    trans_mask = (bin_chr != bin_chr[None].T)
    if mode=='trans':
        mask = trans_mask
    if mode=='cis':
        mask = np.zeros((n,n),dtype=bool)
        if diags is None:
            diags=[1,n]
        for i in range(diags[0],diags[1]+1):
            mask[(range(0,n-i),range(i,n))] = True
            mask[(range(i,n),range(0,n-i))] = True
        mask &= (~trans_mask)
                
    return mask





def gc_strength(A,expected_A,gc,mask=None,d=4):
    '''
    calculate genomic compartment strength of matrix A with respect to gc vector.

    mask: matrix mask used for selecting only a subset of values (see create_mask())

    calculation of gc strength:
    1. calculate matrix LOE by taking log(A/expected_A)
    2. sort gc vector, find 1/d top values and 1/d bottom values
    3. gc_strength = mean(LOE top 1/d vs top 1/d interactions) + mean(LOE bottom 1/d vs bottom 1/d interactions) - 2 * mean(LOE top 1/d vs bottom 1/d interactions)

    '''

    n=A.shape[0]
    LOE = np.log( (A + 1e-7) / (expected_A + 1e-7) )
    LOE = np.clip(LOE,None,np.nanpercentile(LOE,99.9))
    s = np.argsort(gc)

    if mask is not None:
        LOE[~mask]=np.nan

    sorted_LOE=LOE[s,:][:,s]

    nancount=np.sum(np.isnan(gc))
    if nancount>0:
        sorted_LOE = sorted_LOE[:-nancount,:-nancount]
        s = s[:-nancount]

    n=sorted_LOE.shape[0]
    q=int(n/d)

    result = np.nanmean(sorted_LOE[:q,:q]) + np.nanmean(sorted_LOE[-q:,-q:]) - np.nanmean(sorted_LOE[:q,-q:]) - np.nanmean(sorted_LOE[-q:,:q])

    return result, s[:q], s[-q:], LOE, sorted_LOE
