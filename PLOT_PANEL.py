#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 00:53:47 2021

@author: ariad
"""

import pickle, bz2, gzip
from statistics import mean, variance
from math import log
from itertools import starmap

HOME = '/home' #  '/Users' 

def chr_length(chr_id):
    """ Return the chromosome length for a given chromosome, based on the reference genome hg38.""" 
    #The data of chromosome length was taken from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38
    length_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                  'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                  'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                  'chr16': 90338345, 'chr17':  83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                  'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    return length_dict[chr_id]

def std_of_mean(variances):
    """ Standard error of the mean of uncorrelated variables, based on the
        Bienaymé formula. """
    return sum(variances)**.5/len(variances)

def mean_and_var(data):
    """ Calculates the mean and variance. """
    m = mean(data)
    var = variance(data, xbar=m)
    return m, var 

def LLR(y,x):
    """ Calculates the logarithm of y over x and deals with edge cases. """
    if x and y:
        result = log(y/x)
    elif x and not y:
        result = -1.23456789 
    elif not x and y:
        result = +1.23456789 
    elif not x and not y:
        result = 0 
    else:
        result = None    
    return result

def load_likelihoods(filename):
    """ Loads from a file a dictionary that lists genomic windows that contain
    at least two reads and gives the bootstrap distribution of the 
    log-likelihood ratios (LLRs). """
    
    Open = {'bz2': bz2.open, 'gzip': gzip.open}.get(filename.rpartition('.')[-1], open)
    
    with Open(filename, 'rb') as f:
        likelihoods = pickle.load(f)
        info = pickle.load(f)
    return likelihoods, info

def show_info(filename,info,pair):
    S = info['statistics']
    X = info['statistics']['LLRs_per_chromosome'][pair]
    print('\nFilename: %s' % filename)
    print('Depth: %.2f, Chromosome ID: %s, Mean and standard error of meaningful reads per genomic window: %.1f, %.1f.' % (info['depth'], info['chr_id'], S['reads_mean'], S['reads_std']))
    print('Number of genomic windows: %d, Mean and standard error of genomic window size:  %d, %d.' % (S['num_of_windows'],S['window_size_mean'],S['window_size_std']))
    print('Mean LLR: %.3f, Standard error of the mean LLR: %.3f' % ( X['mean'],  X['std_of_mean']))
    print('Fraction of genomic windows with a negative LLR: %.3f' % (X['fraction_of_negative_LLRs']))
    print('The calculation was done in %.3f sec.' % info['runtime'])

def coordinates(windows,chr_id,num_of_buckets):
    """ Lists the buckets and gives the genomic windows that they contain. """
    bin_size = chr_length(chr_id) / num_of_buckets
    result = {}
    j = 0
    
    for i in range(num_of_buckets): ### All buckets before the first the genomic window are filled with Nones.
        if sum(windows[0])/2 < (i+1)*bin_size:
            break
        result[i/num_of_buckets,(i+1)/num_of_buckets] = None
    
    for k,(a,b) in enumerate(windows):
        if not bin_size*i <= (a+b)/2 < bin_size*(i+1):
            result[i/num_of_buckets,(i+1)/num_of_buckets] = (j,k)
            j = k
            for i in range(i+1,num_of_buckets): #Proceed to the next non-empty bucket; Empty buckets are filled with Nones.
                if (a+b)/2 < (i+1)*bin_size:
                    break
                result[i/num_of_buckets,(i+1)/num_of_buckets] = None
    
    for i in range(i,num_of_buckets): ### All buckets after the last the genomic window are filled with Nones.
        result[i/num_of_buckets,(i+1)/num_of_buckets] = (j,k) if j != k else None
        j = k 
    return result

def confidence(LLR_stat,info,num_of_buckets,z_score):
    """ Binning is applied by aggregating the mean LLR of a window across N
        consecutive windows. The boundaries of the bins as well as the mean LLR
        and the standard-error per bin are returned. """
             
    K,M,V = tuple(LLR_stat.keys()), *zip(*LLR_stat.values())
            
    buckets = coordinates(windows=K,chr_id=info['chr_id'],num_of_buckets=num_of_buckets)
    
    X = [*buckets]
    Y = [mean(M[P[0]:P[1]]) if P else None for P in buckets.values()]
    E = [z_score * std_of_mean(V[P[0]:P[1]]) if P else None for P in buckets.values()] 

    return X,Y,E

def detect_transition(X,Y,E):
    A = [(l,j,k) for i,j,k in zip(X,Y,E) for l in i if j!=None and (abs(j)-k)>0 ]
    if len(A)>2:
        x,y,e = zip(*A)
        result = [(x[i+1]+x[i])/2 for i,(a,b) in enumerate(zip(y[:-1],y[1:])) if b/a<0]
    else:
        result = []
        
    return result

def panel_plot(DATA,num_of_buckets_in_chr21,**kwargs):
    scale = 0.5
    z_score = 1
    fs=28 * scale
    import matplotlib as mpl
    save = kwargs.get('save', '')
    if save!='':
            mpl.use('Agg')
    else:
        #['GTK3Agg', 'GTK3Cairo', 'MacOSX', 'nbAgg', 'Qt4Agg', 'Qt4Cairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']
        mpl.use('Qt5Agg')
    mpl.rcParams.update({'figure.max_open_warning': 0})
    import matplotlib.pyplot as plt
    #from scipy.interpolate import interp1d
    num_of_buckets = {'chr'+str(i): num_of_buckets_in_chr21*chr_length('chr'+str(i))//chr_length('chr21') for i in [*range(1,23)]+['X','Y']}

    colors = {frozenset(('BPH','DISOMY')):(177/255,122/255,162/255),
              frozenset(('DISOMY','SPH')):(242/255,142/255,44/255),
              frozenset(('SPH','MONOSOMY')):(239/255,106/255,92/255),
              frozenset(('DISOMY','MONOSOMY')):(104/255,162/255,183/255),
              frozenset(('BPH','SPH')):(104/255,162/255,104/255)}


    fig,axs =plt.subplots(4,6, sharex='col', sharey='row', figsize=(40 * scale, 22.5 * scale))
    #fig.suptitle(kwargs.get('title', ''), fontsize=16)
    #fig.text(0.5, 0, 'Chromosomal position (normalized)',fontsize=28, ha='center')
    #fig.text(0, 0.5, 'Log-likelihood ratio per genomic window', fontsize=28, va='center', rotation='vertical')


    AX = [i for j in axs for i in j]
    
    H = {}
    a,b = 'BPH', 'SPH'
    mean_genomic_window_size = []
    for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA)):

        if (a,b) in info['statistics']['LLRs_per_genomic_window']:
            LLR_stat = info['statistics']['LLRs_per_genomic_window'][(a,b)]
        else:
            _ = {};
            LLR_stat = {window:  mean_and_var([*starmap(LLR, ((_[a], _[b]) for _['MONOSOMY'], _['DISOMY'], _['SPH'], _['BPH'] in L))])
                   for window,L in likelihoods.items()}
        
        X,Y,E = confidence(LLR_stat,info,num_of_buckets=num_of_buckets[info['chr_id']],z_score=1)
        Y = [(y if y else 0) for y in Y]
        E = [(e if e else 0) for e in E]
        mean_genomic_window_size.append(info['statistics']['window_size_mean']/chr_length(info['chr_id']) )
        
        T = [(x[1]+x[0])/2 for x in X]                
        steps_x = [X[0][0]]+[i[1] for i in X[:-1] for j in (1,2)]+[X[-1][1]]
        steps_y = [i for i in Y for j in (1,2)]
        H[a,b] = ax1.plot(steps_x, steps_y, label=f'LLR of {a:s} to {b:s}',color=colors[frozenset((a,b))], linewidth=2, zorder=10, scalex=True, scaley=True, alpha=0.8)
        
        P = [(x[1]-x[0])/2 for x in X]                
        ax1.errorbar(T, Y, xerr = P, ecolor=colors[frozenset((a,b))],marker=None, ls='none',alpha=1, zorder=13, linewidth=5*scale) 
        
        ax1.tick_params(axis='x', labelsize=fs) 
        ax1.tick_params(axis='y', labelsize=fs)
        ax1.xaxis.set_tick_params(width=2*scale)
        ax1.yaxis.set_tick_params(width=2*scale)
        ###ax1.grid(color='black', linestyle='-.', linewidth=1,alpha=0.5)
        for axis in ['top','bottom','left','right']:
            ax1.spines[axis].set_linewidth(2*scale)
        ax1.set_title(info['chr_id'].replace('chr', 'Chromosome '),fontsize=fs)
              
    for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA)):
        ymin,ymax = ax1.get_ylim()
        ym = abs(ymax) if abs(ymax)>abs(ymin) else abs(ymin)        
        ax1.set_ylim((-ym,+ym))
            
    for g,(ax1,(likelihoods,info)) in enumerate(zip(AX,DATA)):
        LLR_stat = info['statistics']['LLRs_per_genomic_window'][(a,b)]
        X0,Y0,E0 = confidence(LLR_stat,info,num_of_buckets=num_of_buckets[info['chr_id']],z_score=z_score)
        Y = [(y if y else 0) for y in Y0]
        E = [(e if e else 0) for e in E0]
        T = [(x[1]+x[0])/2 for x in X0]         
        
        xmin,xmax = ax1.get_xlim()
        ymin,ymax = ax1.get_ylim()
        ax1.errorbar(T, Y, yerr = E, ecolor='black',marker=None, ls='none',alpha=0.2, zorder=15, linewidth=4*scale) 
        ax1.errorbar( xmin + 0.88*(xmax-xmin)-mean_genomic_window_size[g], ymin + 0.12*(ymax-ymin),marker=None, ls='none', xerr=25*mean_genomic_window_size[g], linewidth=2*scale, color='k', capsize=4*scale)
        ax1.text(     xmin + 0.88*(xmax-xmin)-mean_genomic_window_size[g], ymin + 0.09*(ymax-ymin), '25 GW',  horizontalalignment='center', verticalalignment='top',fontsize=2*fs//3)
        ax1.plot([xmin,xmax],[0,0],color='black', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        
        #for i in detect_transition_temp(LLR_stat,chr_length(info['chr_id']),z_score=80):
        for i in detect_transition(X0,Y0,E0):
            #print(g+1,N,len(X0),i)
            ax1.plot([i,i],[ymin,ymax],color='purple', ls='dotted',alpha=0.7,zorder=0, linewidth=2*scale, scalex=False, scaley=False)
        
        ax1.set_xlim((xmin,xmax))                
        ax1.set_ylim((ymin,ymax))
    

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    ###plt.title(f'{RATIO[0]:s} vs. {RATIO[1]:s}\n', fontsize=int(1.2*fs))
    
    plt.xlabel('Chromosomal position (normalized)', fontsize=fs,labelpad=25*scale)
    plt.ylabel('Log-likelihood ratio (normalized)', fontsize=fs,labelpad=45*scale)        
    
    for l in range(1,len(AX)-len(DATA)+1):
        AX[-l].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False, width=0)
        for axis in ['top','bottom','left','right']:
            AX[-l].spines[axis].set_visible(False)
        AX[-l-6].xaxis.set_tick_params(labelbottom=True)
    AX[-1].legend(handles=[i[0] for i in H.values()], title='', bbox_to_anchor=(.5, .45), loc='upper center', ncol=len(H), fancybox=True,fontsize=fs)
    if kwargs.get('title',None): AX[-1].set_title(kwargs.get('title',None), fontsize=int(fs), y=0.55, color='purple')
            
        
    if save!='':
        print('Saving plot...')
        #ax1.set_title(save.rpartition('/')[-1].removesuffix('.png'))
        plt.tight_layout()
        extension = 'svg'
        plt.savefig('.'.join([save,extension]), format=extension, bbox_inches='tight')
        plt.close(fig)
    else:
       plt.tight_layout()
       plt.show() 
       
def main(identifier):      
    #DATA = {f'{identifier:s}.chr{str(i):s}': load_likelihoods(f'/home/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_ZOUVES/{identifier:s}.chr{str(i):s}.LLR.p.bz2') for i in [*range(1,23)]+['X']}
    import matplotlib as mpl
    DATA = {f'{identifier:s}.chr{str(i):s}': load_likelihoods(f'/Users/ariad/Dropbox/postdoc_JHU/LD-PGTA_ecosystem/LD-PGTA_V2/results_RPL/{identifier:s}.chr{str(i):s}.LLR.p.bz2') for i in [*range(1,23)]}
    
    for f,(likelihoods,info) in DATA.items():
        show_info(f'{f:s}.LLR.p.bz2',info,('BPH','SPH'))
    panel_plot(DATA.values(),num_of_buckets_in_chr21=5,title=f'{identifier:s}',save='')

#identifier = '13094FA-B6MTL_3'